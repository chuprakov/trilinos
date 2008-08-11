#include "SundanceMeshReaderBase.hpp"
#include "SundanceExceptions.hpp"
#include "SundancePathUtils.hpp"
#include "SundanceOut.hpp"


using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;


MeshReaderBase::MeshReaderBase(const ParameterList& params)
  : MeshSourceBase(params), 
    filename_()
{
  filename_ = params.get<string>("Filename");
}



int MeshReaderBase::atoi(const std::string& x) const 
{
#ifndef TFLOP
  return std::atoi(x.c_str());
#else
  return ::atoi(x.c_str());
#endif
}

double MeshReaderBase::atof(const std::string& x) const 
{
#ifndef TFLOP
  return std::atof(x.c_str());
#else
  return ::atof(x.c_str());
#endif
}

bool MeshReaderBase::isEmptyLine(const std::string& x) const 
{
  return x.length()==0 || StrUtils::isWhite(x);
}

bool MeshReaderBase::getNextLine(std::istream& is, string& line,
                                         Array<string>& tokens,
                                         char comment) const 
{
  bool rtn = false;
  while (rtn=StrUtils::readLine(is, line))
    {
      SUNDANCE_OUT(this->verbosity() == VerbHigh,
                   "read line [" << line << "]");

      if (line.length() > 0) line = StrUtils::before(line,comment);
      if (isEmptyLine(line)) continue;
      if (line.length() > 0) break;
    }
  tokens = StrUtils::stringTokenizer(line);
  return rtn;
}

RefCountPtr<std::ifstream> MeshReaderBase::openFile(const string& fname, 
                                               const string& description) const
{
  string f = searchForFile(fname);
  RefCountPtr<std::ifstream> rtn = rcp(new std::ifstream(f.c_str()));

  SUNDANCE_OUT(this->verbosity() > VerbMedium,
               "trying to open " << description << " file " << f);

  TEST_FOR_EXCEPTION(rtn.get()==0 || *rtn==0, RuntimeError, 
                     "MeshReaderBase::openFile() unable to open "
                     << description << " file " << f);

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
               "reading " << description << " from " << fname);

  return rtn;
}
