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

#include "Teuchos_ParameterEntryXMLConverter.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"

namespace Teuchos{

RCP<ParameterEntry>
ParameterEntryXMLConverter::fromXMLtoParameterEntry(const XMLObject &xmlObj) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(xmlObj.getRequired(getTypeAttributeName()) != getTypeAttributeValue(),
    BadParameterEntryXMLConverterTypeException,
    "Error: this Parameter Entry XML tag has a type different than "
    "the XMLConverter being used to convert it." <<std::endl <<
    "Parameter name: " << xmlObj.getRequired(
    XMLParameterListWriter::getNameAttributeName()) << std::endl << 
    "XML Parameter Entry type: " << xmlObj.getRequired(getTypeAttributeName()) << std::endl << 
    "XMLConverter type: " << getTypeAttributeValue() << std::endl <<std::endl);
  #endif

  TEST_FOR_EXCEPTION(
    !xmlObj.hasAttribute(getValueAttributeName()), 
    NoValueAttributeExecption,
    ParameterEntry::getTagName() <<" tags must "
    "have a " << getValueAttributeName() << " attribute" << std::endl <<
    "Bad Parameter: " << 
    xmlObj.getAttribute(XMLParameterListWriter::getNameAttributeName()) <<
    std::endl << std::endl);

  RCP<ParameterEntry> toReturn = rcp(new ParameterEntry);
  bool isDefault = false;
  bool isUsed = false;


  if(xmlObj.hasAttribute(getDefaultAttributeName())){
    isDefault = xmlObj.getRequiredBool(getDefaultAttributeName());
  }

  if(xmlObj.hasAttribute(getUsedAttributeName())){
    isUsed = xmlObj.getRequiredBool(getUsedAttributeName());
  }

  setEntryValue(toReturn, xmlObj, isDefault);
  
  if(isUsed){
    toReturn->getAny();
  }
  
  return toReturn;
}


XMLObject
ParameterEntryXMLConverter::fromParameterEntrytoXML(
  RCP<const ParameterEntry> entry, 
  const std::string &name,
  const ParameterEntry::ParameterEntryID& id,
  const XMLParameterListWriter::ValidatorIDsMap& validatorIDsMap) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    (entry->getAny().typeName() != getTypeAttributeValue()) 
    &&
    (
      getTypeAttributeValue() != 
      ParameterEntryXMLConverterDB::getDefaultConverter()->getTypeAttributeValue()
    ),
    BadParameterEntryXMLConverterTypeException,
    "Error: This converter can't convert the given ParameterEntry to XML "
    "because their types don't match." << std::endl <<
    "Parameter name: " << name << std::endl <<
    "Parameter type: " << entry->getAny().typeName() << std::endl <<
    "Converter type: " << getTypeAttributeValue() << std::endl << std::endl);
  #endif

  XMLObject toReturn(ParameterEntry::getTagName());
  toReturn.addAttribute(
    XMLParameterListWriter::getNameAttributeName(), name);
  toReturn.addAttribute(getTypeAttributeName(), getTypeAttributeValue());
  toReturn.addAttribute(getIdAttributeName(), id);
  toReturn.addAttribute(
    getValueAttributeName(), getValueAttributeValue(entry));
  toReturn.addBool(getDefaultAttributeName(), entry->isDefault());
  toReturn.addBool(getUsedAttributeName(), entry->isUsed());
  if(nonnull(entry->validator())){
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(entry->validator()) == validatorIDsMap.end(),
      MissingValidatorDefinitionException,
      "Could not find validator in given ValidatorIDsMap! " << 
      std::endl << std::endl);
    toReturn.addAttribute(
      ValidatorXMLConverter::getIdAttributeName(), 
      validatorIDsMap.find(entry->validator())->second);
  }
  return toReturn;
}


} // namespace Teuchos

