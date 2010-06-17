// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Optika_StandardConditions.hpp"

namespace Optika{

ParameterCondition::ParameterCondition(Teuchos::ParameterEntry* parameter, bool whenParamEqualsValue):
	parameter(parameter),
	whenParamEqualsValue(whenParamEqualsValue){}

BinaryLogicalCondition::BinaryLogicalCondition(std::vector<Teuchos::RCP<Condition> >& conditions):
	conditions(conditions)
{
	if(conditions.size()<2){
		throw InvalidConditionException("Bollocks! Looks like you tried to make a "
		"BinaryLogicalCondition with less than two parameters. You should probably "
		"fix that. \n\n"
		"Error: BinaryLogicalConditions must have at least two conditions");
	}
}

void BinaryLogicalCondition::addCondition(Teuchos::RCP<Condition> toAdd){
	return conditions.push_back(toAdd);
}

OrCondition::OrCondition(std::vector<Teuchos::RCP<Condition> >& conditions):
	BinaryLogicalCondition(conditions){}

bool OrCondition::isConditionTrue(){
	std::vector<Teuchos::RCP<Condition> >::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn || (*it)->isConditionTrue());
	}
	return toReturn;
}

AndCondition::AndCondition(std::vector<Teuchos::RCP<Condition> >& conditions):
	BinaryLogicalCondition(conditions){}

bool AndCondition::isConditionTrue(){
	std::vector<Teuchos::RCP<Condition> >::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn && (*it)->isConditionTrue());
	}
	return toReturn;
}

EqualsCondition::EqualsCondition(std::vector<Teuchos::RCP<Condition> >& conditions):
	BinaryLogicalCondition(conditions){}

bool EqualsCondition::isConditionTrue(){
	std::vector<Teuchos::RCP<Condition> >::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn == (*it)->isConditionTrue());
	}
	return toReturn;
}

NotCondition::NotCondition(Teuchos::RCP<Condition> condition):
	condition(condition){}

bool NotCondition::isConditionTrue(){
	return (!condition->isConditionTrue());
}

StringCondition::StringCondition(Teuchos::ParameterEntry* parameter, std::string value, bool whenParamEqualsValue):
	ParameterCondition(parameter, whenParamEqualsValue), values(ValueList(1,value))
{
	if(!parameter->isType<std::string>()){
		throw InvalidConditionException("The parameter of a String Condition "
		"must be of type string. \n"
		"Expected type: std::string\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

StringCondition::StringCondition(Teuchos::ParameterEntry* parameter, ValueList values, bool whenParamEqualsValue):
	ParameterCondition(parameter, whenParamEqualsValue), values(values)
{
	if(!parameter->isType<std::string>()){
		throw InvalidConditionException("The parameter of a String Condition "
		"must be of type string. \n"
		"Expected type: std::string\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

bool StringCondition::isConditionTrue(){
	bool toReturn = find(values.begin(), values.end(), Teuchos::getValue<std::string>(*parameter)) != values.end();
	return whenParamEqualsValue ? toReturn : !toReturn;
}

BoolCondition::BoolCondition(Teuchos::ParameterEntry* parameter, bool whenParamEqualsValue):
	ParameterCondition(parameter, whenParamEqualsValue)
{
	if(!parameter->isType<bool>()){
		throw InvalidConditionException("The parameter of a Bool Condition "
		"must be of type Bool. \n"
		"Expected type: Bool\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

bool BoolCondition::isConditionTrue(){
	bool toReturn = Teuchos::getValue<bool>(*parameter);
	return whenParamEqualsValue ? toReturn : !toReturn;
}

}

