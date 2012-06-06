// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
/*! \file Zoltan2_Environment.cpp
    \brief The definition of the Environment object.
*/

#ifndef _ZOLTAN2_ENVIRONMENT_CPP_
#define _ZOLTAN2_ENVIRONMENT_CPP_

#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Util.hpp>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>

#include <sstream>
#include <ostream>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////
// Namespace definitions used by this class.

/*! \brief Create an output manager for debugging or status information
 *
 *  \param rank  the MPI rank of the calling process in the application
 *  \param iPrint   true if this process should output information
 *  \param fname    name of file to which output is to be appended, or
 *                      or Z2_UNSET_STRING
 *  \param ost     output stream type
 *  \param mgr     on return, a pointer to the created output manager
 */

void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, int ost,
  Teuchos::RCP<DebugManager> &mgr)
{
  MessageOutputLevel lvl = static_cast<MessageOutputLevel>(level);
  OSType os = static_cast<OSType>(ost);

  bool haveFname = (fname != Z2_UNSET_STRING);
  bool haveStreamName = (os != NUM_OUTPUT_STREAMS);
  bool haveLevel = (lvl != NUM_STATUS_OUTPUT_LEVELS);

  if (!haveFname && !haveStreamName && !haveLevel) {
    // The default if the user didn't specify anything.
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cout, BASIC_STATUS));
    return;
  }

  if (!haveLevel)
    lvl = BASIC_STATUS;

  if (haveFname){
    std::ofstream *dbgFile = new std::ofstream;
    if (iPrint){
      std::string newFname;
      addNumberToFileName(rank, fname, newFname);
      try{
        dbgFile->open(newFname.c_str(), std::ios::out|std::ios::trunc);
      }
      catch(std::exception &e){
        throw std::runtime_error(e.what());
      }
    }
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, *dbgFile, lvl));
    return;
  }

  if (os == COUT_STREAM)
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cout, lvl));
  else if (os == CERR_STREAM)
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cerr, lvl));
  else if (os == NULL_STREAM)
    mgr = Teuchos::rcp(new DebugManager(rank, false, std::cout, lvl));
  else   // should never happen
    throw std::logic_error("invalid debug output stream was not caught");
}

//////////////////////////////////////////////////////////////////////
// Environment definitions

Environment::Environment( Teuchos::ParameterList &problemParams,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
  myRank_(comm->getRank()), numProcs_(comm->getSize()), comm_(comm), 
  errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_(problemParams), params_(problemParams),
  debugOut_(), timerOut_(), timingOn_(false), memoryOut_(), memoryOn_(false)
{
  commitParameters();
}

Environment::Environment():
  myRank_(0), numProcs_(1), comm_(), errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_("emptyList"), params_("emptyList"), 
  debugOut_(), timerOut_(), timingOn_(false), memoryOut_(), memoryOn_(false)
{
  comm_ = Teuchos::DefaultComm<int>::getComm();
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();

  commitParameters();
}

Environment::~Environment()
{
}

void Environment::commitParameters()
{
  using std::string;
  using Teuchos::Array;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  bool emptyList = (params_.begin() == params_.end());

  if (!emptyList){

    ParameterList validParams;
    createValidatorList(params_, validParams);
  
    // Note that since validParams only contains parameters that
    // appear in params_, there are no defaults being set.  We
    // call ParameterList::validateParametersAndSetDefaults() instead of
    // ParameterList::validateParameters() because we want the
    // the validators' validateAndModify() to be called instead
    // of validate().  validateAndModify() "fixes" some of the
    // parameters for us.
  
    params_.validateParametersAndSetDefaults(validParams);

    // For all of the string to integer parameters, convert
    // them to the integer.  I would have
    // expected validateAndModify() to do this.
  
    convertStringToInt(params_);
  }

  /////////////////////////////////////////////////////////////////////

  // Set up for debugging/status output.
  //   By default: if no output stream is specified, then node zero
  //       outputs BASIC_STATUS to std::cout.

#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
  int &level = params_.get<int>("debug_level", NUM_STATUS_OUTPUT_LEVELS);
  string &fname = params_.get<string>("debug_output_file", Z2_UNSET_STRING);
  int &os = params_.get<int>("debug_output_stream", NUM_OUTPUT_STREAMS);

  bool iPrint = (myRank_ == 0);   // default
  const Array<int> *reporters = 
    params_.getPtr<Array<int> >("debug_procs");
  if (reporters)
    iPrint = IsInRangeList(myRank_, *reporters);

  try{
    makeDebugManager(myRank_, iPrint, level, fname, os, debugOut_);
  }
  catch (std::exception &e){
    std::ostringstream oss;
    oss << myRank_ << ": unable to create debug output manager";
    oss << " (" << e.what() << ")";
    throw std::runtime_error(oss.str());
  }
#endif

  // Set up for memory usage output. 
  
#ifndef Z2_OMIT_ALL_PROFILING
  string &f2 = 
    params_.get<string>("memory_output_file", Z2_UNSET_STRING);
  int &os2 = 
    params_.get<int>("memory_output_stream", NUM_OUTPUT_STREAMS);

  const Array<int> *reporters2 = 
    params_.getPtr<Array<int> >("memory_procs");

  bool doMemory = true;

  if (f2 != Z2_UNSET_STRING || os2 != NUM_OUTPUT_STREAMS || reporters2 != NULL){
    // user indicated they want memory usage information
    long numKbytes = 0;
    if (myRank_ == 0)
      numKbytes = getProcessKilobytes();

    Teuchos::broadcast<int, long>(*comm_, 0, 1, &numKbytes);

    if (numKbytes == 0){
      // This is not a Linux system with proc/pid/statm.
      f2 = Z2_UNSET_STRING;
      os2 = NUM_OUTPUT_STREAMS;
      reporters2 = NULL;
      this->debug(BASIC_STATUS, 
        string("Warning: memory profiling requested but not available."));
      doMemory = false;   // can't do it
    }
  }
  else{
    doMemory = false;   // not requested
  }

  if (doMemory){
    bool iPrint = (myRank_ == 0);   // default
    if (reporters2)
      iPrint = IsInRangeList(myRank_, *reporters2);

    try{
      makeMetricOutputManager<long>(myRank_, iPrint, f2, os2, memoryOut_,
        string("KB"), 10);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create memory profiling output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }

    memoryOn_ = true;
  }
#endif

#ifdef Z2_OMIT_ALL_ERROR_CHECKING
  errorCheckLevel_ = NO_ASSERTIONS;
#else
  errorCheckLevel_ = static_cast<AssertionLevel>( 
    params_.get<int>("error_check_level", BASIC_ASSERTION));
#endif
}
  
bool Environment::hasSublist(const Teuchos::ParameterList &pl, 
                             const std::string &listName) const
{
  const Teuchos::ParameterEntry *entry = pl.getEntryPtr(listName);
  if (entry && entry->isList())
    return true;

  return false;
}

void Environment::convertStringToInt(Teuchos::ParameterList &params)
{
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using std::string;
  ParameterList::ConstIterator next = params.begin();

  // Data type of these parameters will now change from string to int

  string validatorName("StringIntegralValidator(int)");
  typedef Teuchos::StringToIntegralParameterEntryValidator<int> s2i_t;

  while (next != params.end()){

    const std::string &name = next->first;
    ParameterEntry &entry = params.getEntry(name);

    if (entry.isList()){
      ParameterList &pl = entry.getValue<ParameterList>(&pl);
      convertStringToInt(pl);
    }
    else{
      if ((entry.validator()).get()){
        if (entry.validator()->getXMLTypeName() == validatorName){
          string &entryValue = entry.getValue<string>(&entryValue);
          RCP<const s2i_t> s2i =
            Teuchos::rcp_dynamic_cast<const s2i_t>(entry.validator(), true);
          int val = s2i->getIntegralValue(entryValue);
          entry.setValue<int>(val);
        }
      }
    }
    ++next;
  }
}

}  //namespace Zoltan2

#endif
