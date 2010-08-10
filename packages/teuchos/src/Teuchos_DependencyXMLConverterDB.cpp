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

#include "Teuchos_DependencyXMLConverterDB.hpp"
#include "Teuchos_StandardDependencyXMLConverters.hpp"
#include "Teuchos_StandardDependencies.hpp"



namespace Teuchos {


#define ADD_TEMPLATED_NUMBER_DEPS(T) \
  ADD_NUMBER_VISUAL_DEP(T) \
  ADD_RANGE_VALIDATOR_DEP(T) \
  ADD_VALIDATOR_ASPECT_DEP(T)

#define ADD_NUMBER_VISUAL_DEP(T) \
  \
  NumberDependency< T > ##T##NumberVisualDependency; \
  masterMap.insert( \
    ConverterPair(##T##NumberVisualDependency.getTypeAttributeValue(), \
      rcp(new NumberVisualDependencyXMLConverter< T >)));

#define ADD_RANGE_VALIDATOR_DEP(T) \
  \
  RangeValidatorDependency< T > ##T##RangeValidatorDependency; \
  masterMap.insert( \
    ConverterPair(##T##RangeValidatorDependency.getTypeAttributeValue(), \
      rcp(new RangeValidatorDependencyXMLConverter< T >)));

#define ADD_VALIDATOR_ASPECT_DEP(T) \
  \
  NumberValidatorAspectDependency< T > \
  ##T##NumberValidatorAspectDependency; \
  \
  masterMap.insert( \
    ConverterPair( \
      ##T##NumberValidatorAspectDependency.getTypeAttributeValue(), \
      rcp(new NumberValidatorAspectDependencyConverter< T >)));


void DependencyXMLConverterDB::addConverter(
  Dependency& dependency,
  RCP<DependencyXMLConverter> converterToAdd){
  getConverterMap().insert(
    ConverterPair(dependency.getTypeAttributeValue(), converterToAdd));
}


RCP<const DependencyXMLConverter>
DependencyXMLConverterDB::getConverter(const Dependency& dependency){
  ConverterMap::const_iterator it = 
    getConverterMap().find(dependency.getTypeAttributeValue());
  TEST_FOR_EXCEPTION(it != getConverterMap().end(),
    CantFindDependencyConverterException,
    "Could not find a DependencyXMLConverter for a dependency"
  )
  return it->second;
}


RCP<const DependencyXMLConverter>
DependencyXMLConverterDB::getConverter(const XMLObject& xmlObject)
{ 
  std::string dependencyType = xmlObject.getRequired(
    DependencyXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(dependencyType);
  TEST_FOR_EXCEPTION(it != getConverterMap().end(),
    CantFindDependencyConverterException,
    "Could not find a DependencyXMLConverter for a dependency"
  )
  return it->second;
}

XMLObject DependencyXMLConverterDB::convertDependency(
  RCP<const Dependency> dependency)
{
  return getConverter(*dependency)->fromDependencytoXML(dependency);
}
 
RCP<ParameterEntryDependency> DependencyXMLConverterDB::convertXML(
  const XMLObject& xmlObject)
{
  return DependencyXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoDependency(xmlObject);
}

DependencyXMLConverterDB::ConverterMap&
DependencyXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  if(masterMap.size() == 0){

    ADD_TEMPLATED_NUMBER_DEPS(int);
    typedef unsigned int uint;
    ADD_TEMPLATED_NUMBER_DEPS(uint);
    typedef short int sint;
    ADD_TEMPLATED_NUMBER_DEPS(sint);
    typedef unsigned short int usint;
    ADD_TEMPLATED_NUMBER_DEPS(usint);
    typedef long int lint;
    ADD_TEMPLATED_NUMBER_DEPS(lint);
    typedef unsigned long int ulint;
    ADD_TEMPLATED_NUMBER_DEPS(ulint);
    ADD_TEMPLATED_NUMBER_DEPS(double);
    ADD_TEMPLATED_NUMBER_DEPS(float);

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    typedef long long int llint;
    ADD_TEMPLATED_NUMBER_DEPS(llint);
    typedef unsigned long long int ullint;
    ADD_TEMPLATED_NUMBER_DEPS(ullint);
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

    StringValidatorDependency stringValidatorDependency;
    masterMap.insert(
      ConverterPair(stringValidatorDependency.getTypeAttributeValue(), 
        rcp(new StringValidatorDependencyXMLConverter)));

    StringVisualDependency stringVisualDependency;
    masterMap.insert(
      ConverterPair(stringVisualDependency.getTypeAttributeValue(), 
        rcp(new StringVisualDependencyXMLConverter)));

    BoolVisualDependency boolVisualDependency;
    masterMap.insert(
      ConverterPair(boolVisualDependency.getTypeAttributeValue(), 
        rcp(new BoolVisualDependencyXMLConverter)));

    BoolValidatorDependency boolValidatorDependency;
    masterMap.insert(
      ConverterPair(boolValidatorDependency.getTypeAttributeValue(), 
        rcp(new BoolValidatorDependencyXMLConverter)));

    NumberArrayLengthDependency numberArrayLengthDependency;
    masterMap.insert(
      ConverterPair(numberArrayLengthDependency.getTypeAttributeValue(),
        rcp(new NumberArrayLengthDependencyConverter)));

  }
  return masterMap;
}


} // namespace Teuchos
