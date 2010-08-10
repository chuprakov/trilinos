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


#ifndef TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDVALIDATORXMLCONVERTERSL_HPP

/*! \file Teuchos_StandardValidatorXMLConverters.hpp
    \brief A collection of standard ValidatorXMLConverters.
*/

#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_DummyObjectGetter.hpp"


namespace Teuchos {


/** \brief Converts StringToIntegralParameterEntryValidators to and from XML.
 */
template<class IntegralType>
class StringToIntegralValidatorXMLConverter : 
  public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;

  /** \brief . */
  XMLObject convertValidator(
    const RCP<const ParameterEntryValidator> validator) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator > 
  getDummyValidator() const{
    return DummyObjectGetter<
    StringToIntegralParameterEntryValidator<IntegralType> >::getDummyObject();
  }
  #endif
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief . */
  static const std::string& getIntegralValueAttributeName() {
    static const std::string integralValueAttributeName_ = "integralValue";
    return integralValueAttributeName_;
  }

  /** \brief . */
  static const std::string& getStringTagName() {
    static const std::string stringTagName_ = "String";
    return stringTagName_;
  }

  /** \brief . */
  static const std::string& getStringValueAttributeName() {
    static const std::string stringValueAttributeName_ = "stringValue";
    return stringValueAttributeName_;
  }

  /** \brief . */
  static const std::string& getStringDocAttributeName() {
    static const std::string stringDocAttributeName_ = "stringDoc";
    return stringDocAttributeName_;
  }

  /** \brief . */
  static const std::string& getDefaultParameterAttributeName() {
    static const std::string defaultParameterAttributeName_ =
      "defaultParameterName";
    return defaultParameterAttributeName_;
  }
  
  //@}

};


//
// Implementations
//


template<class IntegralType>
RCP<ParameterEntryValidator>
StringToIntegralValidatorXMLConverter<IntegralType>::convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;
{
  Array<std::string> strings;
  Array<std::string> stringDocs;
  Array<IntegralType> integralValues;
  for (int i=0; i<xmlObj.numChildren(); ++i) {
    XMLObject currentChild = xmlObj.getChild(i);
    TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
      BadTagException,  
      "Error converting xmlObject to "
      "StringToIntegralParameterEntryValidator." << std::endl << 
      "Unrecognized tag: " << currentChild.getTag());
    strings.append(currentChild.getRequired(getStringValueAttributeName()));
    if (currentChild.hasAttribute(getIntegralValueAttributeName())) {
      integralValues.append(
        currentChild.getRequired<IntegralType>(
          getIntegralValueAttributeName()));
    }
    if (currentChild.hasAttribute(getStringDocAttributeName())) {
      stringDocs.append(
        currentChild.getRequired<std::string>(getStringDocAttributeName()));
    }
  }
  std::string defaultParameterName = 
    xmlObj.getRequired(getDefaultParameterAttributeName());

  return stringToIntegralParameterEntryValidator<IntegralType>(
    strings, stringDocs, integralValues, defaultParameterName, validatorID);
}


template<class IntegralType>
XMLObject StringToIntegralValidatorXMLConverter<IntegralType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator) const
{
  RCP<const StringToIntegralParameterEntryValidator<IntegralType> > 
    castedValidator =
    rcp_dynamic_cast<
      const StringToIntegralParameterEntryValidator<IntegralType> >(
        validator, true);

  XMLObject toReturn(validator->getXMLTagName());

  RCP<const Array<std::string> > stringValues =
    castedValidator->validStringValues();
  RCP<const Array<std::string> > stringDocValues = 
    castedValidator->getStringDocs();

  bool hasStringDocs = 
    !(stringDocValues.is_null()) && (stringDocValues->size() != 0);
  for (int i =0; i<stringValues->size(); ++i) {
    XMLObject stringTag(getStringTagName());
    stringTag.addAttribute(getStringValueAttributeName(), (*stringValues)[i]);
    stringTag.addAttribute(getIntegralValueAttributeName(),
      castedValidator->getIntegralValue((*stringValues)[i]));
    if (hasStringDocs) {
      stringTag.addAttribute(
        getStringDocAttributeName(), (*stringDocValues)[i]);
    }
    toReturn.addChild(stringTag);
  }
  toReturn.addAttribute(getDefaultParameterAttributeName(),
    castedValidator->getDefaultParameterName());
  toReturn.addAttribute(getIntegralValueAttributeName(),
    TypeNameTraits<IntegralType>::name());
  return toReturn;
}

/*
 * \brief Converts AnyNumberParameterEntryValidators to and from XML.
 */
class AnyNumberValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;

  /** \brief . */
  XMLObject convertValidator(
   const RCP<const ParameterEntryValidator> validator) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief . */
  static const std::string& getAllowIntAttributeName() {
    static const std::string allowIntAttributeName_ = "allowInt";
    return allowIntAttributeName_;
  }

  /** \brief . */
  static const std::string& getAllowDoubleAttributeName() {
    static const std::string allowDoubleAttributeName_ = "allowDouble";
    return allowDoubleAttributeName_;
  }

  /** \brief . */
  static const std::string& getAllowStringAttributeName() {
    static const std::string allowStringAttributeName_ = "allowString";
    return allowStringAttributeName_;
  }
  
  /** \brief . */
  static const std::string& getPrefferedTypeAttributeName() {
    static const std::string prefferedTypeAttributeName_ = "prefferedType";
    return prefferedTypeAttributeName_;
  }
  
  //@}
 
};


/** \brief Converts EnhancedNumberValidators to and from XML.
 */
template<class T>
class EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;

  /** \brief . */
  XMLObject convertValidator(
    const RCP<const ParameterEntryValidator> validator)const;

#ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<EnhancedNumberValidator<T> >::getDummyObject();
  }
#endif
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief . */
  static const std::string& getMinAttributeName() {
    static const std::string minAttributeName = "min";
    return minAttributeName;
  }

  /** \brief . */
  static const std::string& getMaxAttributeName() {
    static const std::string maxAttributeName = "max";
    return maxAttributeName;
  }

  /** \brief . */
  static const std::string& getStepAttributeName() {
    static const std::string stepAttributeName = "step";
    return stepAttributeName;
  }

  /** \brief . */
  static const std::string& getPrecisionAttributeName() {
    static const std::string precisionAttributeName = "precision";
    return precisionAttributeName;
  }
  
  //@}

};


template<class T>
RCP<ParameterEntryValidator> EnhancedNumberValidatorXMLConverter<T>::convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;
{
  RCP<EnhancedNumberValidator<T> > toReturn = 
    rcp(new EnhancedNumberValidator<T>(validatorID));
  xmlObj.getWithDefault(
    getStepAttributeName(), EnhancedNumberTraits<T>::defaultStep()),
  xmlObj.getWithDefault(getPrecisionAttributeName(),
    EnhancedNumberTraits<T>::defaultPrecision());
  if (xmlObj.hasAttribute(getMinAttributeName())) {
    toReturn->setMin(xmlObj.getRequired<T>(getMinAttributeName()));
  }
  if (xmlObj.hasAttribute(getMaxAttributeName())) {
    toReturn->setMax(xmlObj.getRequired<T>(getMaxAttributeName()));
  }
  return toReturn;
}


template<class T>
XMLObject EnhancedNumberValidatorXMLConverter<T>::convertValidator(
  const RCP<const ParameterEntryValidator > validator) const
{
  RCP<const EnhancedNumberValidator<T> > castedValidator =
    rcp_dynamic_cast<const EnhancedNumberValidator<T> >(validator, true);
  XMLObject toReturn(castedValidator->getXMLTagName());
  if (castedValidator->hasMin()) {
    toReturn.addAttribute<T>(getMinAttributeName(), castedValidator->getMin());
  }
  if (castedValidator->hasMax()) {
    toReturn.addAttribute<T>(getMaxAttributeName(), castedValidator->getMax());
  }
  toReturn.addAttribute<T>(getStepAttributeName(), castedValidator->getStep());
  toReturn.addAttribute<T>(
    getPrecisionAttributeName(), castedValidator->getPrecision());
  return toReturn;
}


/*
 * \brief Converts FileNameValidators to and from XML.
 */
class FileNameValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;

  /** \brief . */
  XMLObject convertValidator(
    const RCP<const ParameterEntryValidator> validator) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief . */
  static const std::string& getFileMustExistAttributeName() {
    static const std::string fileMustExistAttributeName = "fileMustExist";
    return fileMustExistAttributeName;
  }
  
  //@}
  
};


/*
 * \brief Converts StringValidators to and from XML.
 */
class StringValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;

  /** \brief . */
  XMLObject convertValidator(
    const RCP<const ParameterEntryValidator> validator) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif
  
  //@}

private:
  
  /** \name Private Members */
  //@{
  
  /** \brief . */
  static const std::string& getStringTagName() {
    static const std::string stringTagName = "String";
    return stringTagName;
  }

  /** \brief . */
  static const std::string& getStringValueAttributeName() {
    static const std::string stringValueAttributeName = "stringValue";
    return stringValueAttributeName;
  }
  
  //@}

};


/*
 * \brief Converts ArrayValidators to and from XML.
 */
template<class ValidatorType, class EntryType>
class ArrayValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;

  /** \brief . */
  XMLObject convertValidator(
    const RCP<const ParameterEntryValidator> validator) const;

#ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<ArrayValidator<ValidatorType, EntryType> >::
      getDummyObject();
  }
#endif
  
  //@}

};


template<class ValidatorType, class EntryType>
RCP<ParameterEntryValidator>
ArrayValidatorXMLConverter<ValidatorType, EntryType>::convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;
{
  RCP<ValidatorType> prototypeValidator;
  if (xmlObj.hasAttribute(getPrototypeIdAttributeName())) {
    RCP<ParameterEntryValidator> result =
      ParameterEntryValidator::getValidator(
        xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
          getPrototypeIdAttributeName()));
    if (result != validatorMap.end()) {
      prototypeValidator = 
        rcp_dynamic_cast<ValidatorType>(result->second, true);
    }
    else {
      TEST_FOR_EXCEPTION(true,
        std::runtime_error,
        "Could not find prototype validator with id: "
        << xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
          getPrototypeIdAttributeName()) << std::endl<< std::endl);
    }
  }
  else {
    prototypeValidator = rcp_dynamic_cast<ValidatorType>(
      ValidatorXMLConverterDB::convertXML(
        xmlObj.getChild(0), validatorMap), true);
  }
  return rcp(new ArrayValidator<ValidatorType, EntryType>(prototypeValidator));
}


template<class ValidatorType, class EntryType>
XMLObject 
ArrayValidatorXMLConverter<ValidatorType, EntryType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator) const
{
  XMLObject toReturn(validator->getXMLTagName());
  RCP<const ArrayValidator<ValidatorType, EntryType> > castedValidator = 
    rcp_dynamic_cast<const ArrayValidator<ValidatorType, EntryType> >(
      validator, true);
  ParameterEntryValidator::ValidatorID prototypeID = 
    ParameterEntryValidator::getValidatorID(castedValidator->getPrototype()) 
  if (
    prototypeID
    != 
    OrdinalTraits<ParameterEntryValidator::ValidatorID>::invalid())
  {
    toReturn.add<ParameterEntryValidator::ValidatorID>(
      getPrototypeIdAttributeName(), prototypeID);
  }
  else {
    toReturn.addChild(
    ValidatorXMLConverterDB::convertValidator(castedValidator->getPrototype()));
  }
  return toReturn;
}

} // namespace Teuchos


#endif  // TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP

