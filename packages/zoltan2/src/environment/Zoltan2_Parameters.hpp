// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_PARAMETERS_HPP_
#define _ZOLTAN2_PARAMETERS_HPP_

#include <Zoltan2_config.h>

#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cctype>

#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>

/*! \file Zoltan2_Parameters.hpp

  This file contains parameter-related definitions.
*/

// Had to redefine this type from Teuchos_ParameterEntryValidator.hpp.
// Compiler stumbled on it.
typedef Teuchos::RCP<const Teuchos::Array<std::string> > ValidStringsList;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////

template <typename Integral>
  class IntegerRangeListValidator : public Teuchos::ParameterEntryValidator
{
private:
  Integral min_;
  Integral max_;

  static const std::string listDelim_;
  static const std::string rangeDelim_;
  static const std::string allText_;

  static void checkValid(char c); 
  static bool listSaysAll(std::string &l);
  static int breakRange(std::string &range, std::string &from, std::string &to);

public:
  // Constructor: any Integral is valid
  IntegerRangeListValidator();

  // Constructor: only Integrals in [validMin,validMax] are valid
  IntegerRangeListValidator(Integral validMin, Integral validMax); 

  // Implementation of ParameterEntryValidator interface

  const std::string getXMLTypeName() const; 

  void printDoc(std::string const& docString, std::ostream &out) const;

  //Teuchos::ValidStringsList validStringValues() const ;
  ValidStringsList validStringValues() const ;

  void validate( Teuchos::ParameterEntry  const& entry,
    std::string const& paramName, std::string const& sublistName
    ) const;

  void validateAndModify( std::string const& paramName,
    std::string const& sublistName, Teuchos::ParameterEntry * entry
    ) const;
}; // end class

// Helpers for IntegralRangeList parameter type

template <typename Integral>
  bool validIntegralRangeList(const Teuchos::Array<Integral> &vals);

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::Array<Integral> &vals);

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::ParameterEntry &e);

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::Array<Integral> &vals);

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::ParameterEntry &e);

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e);

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::Array<Integral> &vals);

template <typename Integral>
  void printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl);

// Function declarations

void createValidParameterList(Teuchos::ParameterList &pl, 
  const Teuchos::Comm<int> &comm);

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

// Parameter enumerated types.
//
//  If you change these enumerators, change their documentation
//  in Zoltan2_Parameters.cpp.
//

enum RangeType {
  RANGE_INCLUDES_ALL,
  RANGE_IS_EMPTY,
  RANGE_IS_LISTED,
  NUM_RANGE_TYPES};

enum AssertionLevel {
  BASIC_ASSERTION,    /*!< checks that should always be done (user input) */
  COMPLEX_ASSERTION,  /*!< checks that take extra time (valid input graph?) */
  DEBUG_MODE_ASSERTION,   /*!< only done when check absolutely everything */
  NUM_ASSERTION_LEVELS};

enum MessageOutputLevel {
  NO_STATUS,                 /*!< don't display status/debug messages */
  BASIC_STATUS,              /*!< the status at each high level step */
  DETAILED_STATUS,           /*!< include status at algorithm sub-steps */
  VERBOSE_DETAILED_STATUS,   /*!< include more detail about sub-steps */
  NUM_STATUS_OUTPUT_LEVELS};

enum MessageSummaryLevel{
  LOCAL_SUMMARY,              /*!< messages should display local info only */
  GLOBAL_SUMMARY,             /*!< include global min, max, avg, etc. */
  NUM_STATUS_SUMMARY_LEVELS};

template <typename Integral>
const std::string IntegerRangeListValidator<Integral>::listDelim_(",");

template <typename Integral>
const std::string IntegerRangeListValidator<Integral>::rangeDelim_("-");

template <typename Integral>
const std::string IntegerRangeListValidator<Integral>::allText_("all");

template <typename Integral>
  void IntegerRangeListValidator<Integral>::checkValid(char c)
{
  if (std::isspace(c) || std::isdigit(c) || (c == ',') || (c == '-'))
    return;
  else
    throw std::runtime_error("invalid integer range list");
}

template <typename Integral>
  bool IntegerRangeListValidator<Integral>::listSaysAll(std::string &l)
{
  std::transform(l.begin(), l.end(), l.begin(), tolower);
  if (l.find(allText_) != std::string::npos)
    return true;  // "all" is in the string
  else
    return false;
}

template <typename Integral>
  int IntegerRangeListValidator<Integral>::breakRange(
    std::string &range, std::string &from, std::string &to)
{
  from.clear();
  to.clear();
  std::string::size_type loc = range.find(rangeDelim_);
  if (loc == std::string::npos){
    from = range;
  }
  else{
    from = range.substr(0, loc);
    to = range.substr(loc+1, range.size());
  }
  long a, b;
  std::istringstream iss1(from);
  iss1 >> a;
  b = a;
  if (to.size() > 0){
    std::istringstream iss2(to);
    iss2 >> b;
    if (b < a)
      std::swap(from,to);
  }
  return (b != a) ? 2 : 1;
}


template <typename Integral>
  IntegerRangeListValidator<Integral>::IntegerRangeListValidator(): 
    min_(1), max_(0)
{
}

template <typename Integral>
  IntegerRangeListValidator<Integral>::IntegerRangeListValidator(
    Integral validMin, Integral validMax) :
      min_(validMin), max_(validMax)
{
  if (min_ < max_) std::swap(min_,max_);
}

  // Implementation of ParameterEntryValidator interface

template <typename Integral>
  const std::string 
    IntegerRangeListValidator<Integral>::getXMLTypeName() const 
{
  std::string className("IntegerRangeListValidator");
  std::string classType("("+Teuchos::TypeNameTraits<Integral>::name()+")");
  return std::string(className + classType);
}

template <typename Integral>
  void IntegerRangeListValidator<Integral>::printDoc(
    std::string const& docString, std::ostream &out) const  
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#\tAn integer range list is a string which can contain:\n";
  out << "#\t\tthe text \"all\", which indicates all values\n";
  out << "#\t\ta list of integer ranges separated by commas.\n";
  out << "#\tA range is one value, or two values separated by a dash.\n";
  out << "#\tExample: \"all\" or \"1-10\" or \"3, 10-12\" or \"25\"\n";
  if (max_ >= min_){
    out << "#\tThe range of valid integers is [";
    out << min_ << "," << max_ << "]\n";
  }
}

template <typename Integral>
  ValidStringsList IntegerRangeListValidator<Integral>::validStringValues() const 
{ 
  return Teuchos::null; 
}

template <typename Integral>
  void IntegerRangeListValidator<Integral>::validate( 
    Teuchos::ParameterEntry  const& entry,
    std::string const& paramName, std::string const& sublistName) const
{
  if (!entry.isType<std::string>()){
    return;  // already converted to an an array
  }
  std::string *sptr=NULL;
  std::string &rangeList = entry.getValue(sptr);
  std::string inValue(rangeList);

  if (listSaysAll(inValue))
    return;  // "all" is in the string

  // throw error if invalid integer range list
  std::for_each(inValue.begin(), inValue.end(), checkValid);

  if (max_ >= min_){
    std::string::const_iterator rangeBegin = inValue.begin();
    std::string::const_iterator valueEnd = inValue.end();

    while (rangeBegin != valueEnd){
      std::string::const_iterator rangeEnd = std::search(
        rangeBegin, valueEnd, listDelim_.begin(), listDelim_.end());
      std::string range(rangeBegin, rangeEnd);
      std::string aHalf, bHalf;
      int count = breakRange(range, aHalf, bHalf);

      Integral a, b;
      std::istringstream iss1(aHalf);
      iss1 >> a;
      if (count > 1){
        std::istringstream iss2(bHalf);
        iss2 >> b;
      }
      else
        b = a;

      if ((a < min_) || (b > max_)){
        std::ostringstream oss;
        oss << "input range [" << a << "," << b << "] ";
        oss << "exceeds valid range [" << min_ << "," << max_ << "] ";
        throw std::runtime_error(oss.str());
      }
      if (rangeEnd == valueEnd)
        rangeBegin = rangeEnd;
      else
        rangeBegin = ++rangeEnd;
    }
  }
}

template <typename Integral>
  void IntegerRangeListValidator<Integral>::validateAndModify( 
    std::string const& paramName, std::string const& sublistName, 
    Teuchos::ParameterEntry * entry) const
{
  typedef typename Teuchos::Array<Integral>::size_type arraySize_t;
  if (!entry->isType<std::string>()){
    return;
  }
  
  std::string *sptr=NULL;
  std::string &rangeList = entry->getValue(sptr);
  Teuchos::Array<Integral> valueList;

  std::string inValue(rangeList);

  if (listSaysAll(inValue)){
    valueList.push_back(RANGE_INCLUDES_ALL);
  }
  else{
    // throw error if invalid integer range list
    std::for_each(inValue.begin(), inValue.end(), checkValid);

    std::string::const_iterator rangeBegin = inValue.begin();
    std::string::const_iterator valueEnd = inValue.end();

    while (rangeBegin != valueEnd){
      std::string::const_iterator rangeEnd = std::search(rangeBegin,
        valueEnd, listDelim_.begin(), listDelim_.end());
      std::string range(rangeBegin, rangeEnd);
      std::string aHalf, bHalf;
      int count = breakRange(range, aHalf, bHalf);

      Integral a, b;
      std::istringstream iss1(aHalf);
      iss1 >> a;
      if (count > 1){
        std::istringstream iss2(bHalf);
        iss2 >> b;
      }
      else
        b = a;

      if ((max_ >= min_) && ((a < min_) || (b > max_))){
        std::ostringstream oss;
        oss << "input range [" << a << "," << b << "] ";
        oss << "exceeds valid range [" << min_ << "," << max_ << "] ";
        throw std::runtime_error(oss.str());
      }
      for (Integral i=a; i <=b; i++)
        valueList.push_back(i);
      if (rangeEnd == valueEnd)
       rangeBegin = rangeEnd;
      else
        rangeBegin = ++rangeEnd;
    }
    if (valueList.size() > 1){  // sort & remove duplicates
      std::sort(valueList.begin(), valueList.end());
      arraySize_t listEnd = 0;
      arraySize_t length = valueList.size();
      for (arraySize_t i=1; i < length; i++){
        if (valueList[i] > valueList[listEnd]){
          listEnd++;
          if (listEnd != i)
            valueList[listEnd] = valueList[i];
        }
      }
      if (++listEnd < length)
        valueList.resize(listEnd);
    }

    Integral flag = RANGE_IS_LISTED;
    if (valueList.size() == 0){
      flag = RANGE_IS_EMPTY;
    }
    else if (max_ >= min_){
      Integral allSize = max_ - min_ + 1;
      if (valueList.size() == allSize){
        flag = RANGE_INCLUDES_ALL;
        valueList.clear();
      }
    }
    valueList.push_back(flag);
  }
  entry->setValue(valueList);
}

// Helpers for IntegralRangeList parameter type

template <typename Integral>
  bool validIntegralRangeList(const Teuchos::Array<Integral> &vals)
{
  typedef typename Teuchos::Array<Integral>::size_type arraySize_t;
  arraySize_t len = vals.size();
  if (len==0)
    return false;

  Integral flag = vals[len-1];
  if ((flag != RANGE_INCLUDES_ALL) && (flag != RANGE_IS_EMPTY) &&
      (flag != RANGE_IS_LISTED))
    return false;

  return true;
}

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
{
  if (!validIntegralRangeList(vals))
    throw std::runtime_error("list is not a valid range list");
  return vals[vals.size()-1] == RANGE_INCLUDES_ALL;
}

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &vals = e.getValue(valPtr);
  return allValuesAreInRangeList(vals);
}

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
{
  if (!validIntegralRangeList(vals))
    throw std::runtime_error("list is not a valid range list");
  return vals[vals.size()-1] == RANGE_IS_EMPTY;
}

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &vals = e.getValue(valPtr);
  return noValuesAreInRangeList(vals);
}

// TODO :
// inList(std::vector<Integral> &val, std::vector<bool> &result)

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::Array<Integral> &valList)
{
  if (allValuesAreInRangeList(valList))
    return true;
  else if (noValuesAreInRangeList(valList))
    return false;

  typename Teuchos::Array<Integral>::const_iterator flag = valList.end();
  --flag;
  if (std::binary_search(valList.begin(), flag, val))
    return true;
  else
    return false;
}

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &valList = e.getValue(valPtr);

  return IsInRangeList(val, valList);
} 

template <typename Integral>
  void printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl)
{
  if (Zoltan2::allValuesAreInRangeList(irl))
    os << "all";
  else if (Zoltan2::noValuesAreInRangeList(irl))
    os << "empty";
  else{
    Teuchos::ArrayView<const Integral> view = 
      irl.view(0, irl.size()-1); // skip last value, it's a flag
    os << view;
  }    
}

}  // end of namespace Zoltan2

#endif
