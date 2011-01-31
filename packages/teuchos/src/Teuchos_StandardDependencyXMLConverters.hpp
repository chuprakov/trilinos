// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardDependencyXMLConverters.hpp
 * \brief A collection of standard DependencyXMLConverters.
*/


#include "Teuchos_DependencyXMLConverter.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_XMLDependencyExceptions.hpp"


namespace Teuchos {


/** \brief An xml converter for VisualDepenencies
 */
class VisualDependencyXMLConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific visual dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   * @return An XMLObject representing the VisualDependency.
   */
  virtual void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const = 0;

  /** \brief Converts any special aspects of a
   * specific visual dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependees The dependees of the visual dependency.
   * @param dependents The dependents of the visua dependency.
   * @param showIf The showIf attribute of the visual dependency.
   * @param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   * @return The converted VisualDependency.
   */
  virtual RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  /** \brief Gets the name of the showif attribute */
  static const std::string& getShowIfAttributeName(){
    static const std::string showIfAttributeName = "showIf";
    return showIfAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for ValidatorDependencies.
 */
class ValidatorDependencyXMLConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific validator dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @param validatorIDsMap A map containing validators and their associated
   * IDs.
   * @return An XMLObject representing the ValidatorDepenency.
   */
  virtual void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const = 0;

  /** \brief Converts any special aspects of a
   * specific validator dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependee The dependees of the validator dependency.
   * @param dependents The dependents of the validator dependency.
   * @param validatorIDsMap A map containing validators and their associated
   * IDs.
   * @return The converted ValidatorDependency.
   */
  virtual RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  //@}
  
};

/** \brief An xml converter for StringVisualDepenencies
 */
class StringVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  /** \brief Gets the StringValues Tag */
  static const std::string& getStringValuesTagName(){
    static const std::string stringValuesTagName = "StringValues";
    return stringValuesTagName;
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief Gets the String Tag */
  static const std::string& getStringTagName(){
    static const std::string stringTagName = "string";
    return stringTagName;
  }

  /** \brief Gets the Value attribute name */
  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for BoolVisualDepenencies
 */
class BoolVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  //@}

};

/** \brief An xml converter for NumberVisualDependencies
 */
template<class T>
class NumberVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;
  
  //@}

};

template<class T>
void NumberVisualDependencyXMLConverter<T>::convertSpecialVisualAttributes(
  RCP<const VisualDependency> dependency,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  
}
  
template<class T>
RCP<VisualDependency> 
NumberVisualDependencyXMLConverter<T>::convertSpecialVisualAttributes(
  const XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A NumberVisualDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  return rcp(new NumberVisualDependency<T>(
    *(dependees.begin()), dependents));
}

/** \brief An xml converter for ConditionVisualDependencies
 */
class ConditionVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;
  
  //@}

};


/** \brief An xml converter for StringValidatorDependencies
 */
class StringValidatorDependencyXMLConverter : 
  public ValidatorDependencyXMLConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  static const std::string& getValuesAndValidatorsTag(){
    static const std::string valuesAndValidatorsTag = "ValuesAndValidators";
    return valuesAndValidatorsTag;
  }
  
  //@}

private:
  /** \name Private Members */
  //@{
 
  /** \brief . */
  static const std::string& getPairTag(){
    static const std::string pairTag = "Pair";
    return pairTag;
  }

  /** \brief . */
  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  /** \brief . */
  static const std::string& getValidatorIDAttributeName(){
    static const std::string validatorIDAttributeName = "validatorID";
    return validatorIDAttributeName;
  }

  /** \brief . */
  static const std::string& getDefaultValidatorIDAttributeName(){
    static const std::string defaultValidatorIDAttributeName = 
      "defaultValidatorID";
    return defaultValidatorIDAttributeName;
  }
  
  //@}
  
};

/** \brief An xml converter for BoolValidatorDependencies
 */
class BoolValidatorDependencyXMLConverter : public ValidatorDependencyXMLConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const;
  
  //@}

private:

  /** \name Private Members */
  //@{
    
    /** \brief . */
    static const std::string& getFalseValidatorIDAttributeName(){
      static const std::string falseValidatorIDAttributeName = 
        "falseValidatorID";
      return falseValidatorIDAttributeName;
    }

    /** \brief . */
    static const std::string& getTrueValidatorIDAttributeName(){
      static const std::string trueValidatorIDAttributeName = 
        "trueValidatorID";
      return trueValidatorIDAttributeName;
    }
  
  //@}
  
};

/** \brief An xml converter for RangeValidatorDependencies
 */
template<class T>
class RangeValidatorDependencyXMLConverter : 
  public ValidatorDependencyXMLConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  static const std::string& getRangesAndValidatorsTag(){
    static const std::string rangesAndValidatorsTag = "RangesAndValidators";
    return rangesAndValidatorsTag;
  }
  
  //@}

private:
  /** \name Private Members */
  //@{
 

  /** \brief . */
  static const std::string& getPairTag(){
    static const std::string pairTag = "Pair";
    return pairTag;
  }

  /** \brief . */
  static const std::string& getMinAttributeName(){
    static const std::string minAttributeName = "min";
    return minAttributeName;
  }

  /** \brief . */
  static const std::string& getMaxAttributeName(){
    static const std::string maxAttributeName = "max";
    return maxAttributeName;
  }


  /** \brief . */
  static const std::string& getValidatorIDAttributeName(){
    static const std::string validatorIDAttributeName = "validatorID";
    return validatorIDAttributeName;
  }

  //@}
  
};

template<class T>
void
RangeValidatorDependencyXMLConverter<T>::convertSpecialValidatorAttributes(
  RCP<const ValidatorDependency> dependency,
  XMLObject& xmlObj,
  ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const RangeValidatorDependency<T> > castedDependency = 
    rcp_dynamic_cast<const RangeValidatorDependency<T> >(dependency, true);

  XMLObject rangesAndValidatorsTag(getRangesAndValidatorsTag());

  castedDependency->getRangeToValidatorMap();
  for(
    typename RangeValidatorDependency<T>::RangeToValidatorMap::const_iterator 
    it = castedDependency->getRangeToValidatorMap().begin();
    it != castedDependency->getRangeToValidatorMap().end();
    ++it)
  {
    T min = it->first.first;
    T max = it->first.second;
    if(validatorIDsMap.find(it->second) == validatorIDsMap.end()){
      validatorIDsMap.insert(it->second);
    }
    ParameterEntryValidator::ValidatorID validatorID = 
      validatorIDsMap.find(it->second)->second;
    XMLObject pairTag(getPairTag());
    pairTag.addAttribute(getMinAttributeName(), min);
    pairTag.addAttribute(getMaxAttributeName(), max);
    pairTag.addAttribute(getValidatorIDAttributeName(), validatorID);
    rangesAndValidatorsTag.addChild(pairTag);
  }
  xmlObj.addChild(rangesAndValidatorsTag);
}

template<class T>
RCP<ValidatorDependency> 
RangeValidatorDependencyXMLConverter<T>::convertSpecialValidatorAttributes(
  const XMLObject& xmlObj,
  RCP<const ParameterEntry> dependee,
  const Dependency::ParameterEntryList dependents,
  const IDtoValidatorMap& validatorIDsMap) const
{

  int result = xmlObj.findFirstChild(getRangesAndValidatorsTag()); 
  TEST_FOR_EXCEPTION(result == -1,
    MissingRangesAndValidatorsTagException,
    "Error: All RangeValidatorDependencies must have a " << 
    getRangesAndValidatorsTag() << " tag!" << std::endl << std::endl);

  XMLObject rangesAndValidatorsTag = xmlObj.getChild(result);

  typename RangeValidatorDependency<T>::RangeToValidatorMap 
    rangesAndValidators;
  for(int i = 0 ; i < rangesAndValidatorsTag.numChildren(); ++i){
    XMLObject child = rangesAndValidatorsTag.getChild(i);
    T min = child.getRequired<T>(getMinAttributeName());
    T max = child.getRequired<T>(getMaxAttributeName());
    ParameterEntryValidator::ValidatorID currentID =
      child.getRequired<ParameterEntryValidator::ValidatorID>(
          getValidatorIDAttributeName());
      
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(currentID) == validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find validator in given ValidatorIDsMap! " << std::endl <<
      std::endl);
    RCP<ParameterEntryValidator> validator = 
      validatorIDsMap.find(currentID)->second;
   
    rangesAndValidators.insert(
      typename RangeValidatorDependency<T>::RangeValidatorPair(
        typename RangeValidatorDependency<T>::Range(min, max), validator));
  }

  return rcp(new RangeValidatorDependency<T>(
    dependee, dependents, rangesAndValidators));
}

/** \brief An xml converter for NumberArrayLengthDependencies.
 */
template<class DependeeType, class DependentType>
class NumberArrayLengthDependencyXMLConverter : public DependencyXMLConverter{

public:

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  //@}

};

template<class DependeeType, class DependentType>
RCP<Dependency> 
NumberArrayLengthDependencyXMLConverter<DependeeType, DependentType>::convertXML(
  const XMLObject& xmlObj, 
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap,
  const IDtoValidatorMap& validatorIDsMap) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A NumberArrayLengthDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  
  return rcp(
    new NumberArrayLengthDependency<DependeeType, DependentType>(
      *(dependees.begin()), dependents));
}

template<class DependeeType, class DependentType>
void
NumberArrayLengthDependencyXMLConverter<DependeeType, DependentType>::convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const
{
  
}


} // namespace Teuchos


#endif // TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

