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

ParameterCondition::ParameterCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue):
	Condition(Condition::ParamCon),
	parameterName(parameterName),
	parentList(parentList),
	whenParamEqualsValue(whenParamEqualsValue)
{
	parameter = parentList->getEntryPtr(parameterName);
	if(parameter == NULL){
		throw InvalidConditionException("Oh noes!!!!!!! Looks like the parameter " +
		parameterName + " isn't actually contained in the " + parentList->name() + " parameter list. "
		"You should go back and check your code. Maybe the information below can help you.\n\n"
		"Error: Parameter not contained in specified Parent list:\n"
		"Problem Parameter: " + parameterName + "\n"
		"Problem List: " + parentList->name());
	}
}

bool ParameterCondition::containsAtLeasetOneParameter(){
	return true;
}

Dependency::ParameterParentMap ParameterCondition::getAllParameters(){
	Dependency::ParameterParentMap toReturn;
	toReturn.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(parameterName, parentList));
	return toReturn;
}

BinaryLogicalCondition::BinaryLogicalCondition(ConditionList& conditions):
	Condition(Condition::BinLogicCon),
	conditions(conditions){}

void BinaryLogicalCondition::addCondition(Teuchos::RCP<Condition> toAdd){
	conditions.append(toAdd);
}

bool BinaryLogicalCondition::containsAtLeasetOneParameter(){
	for(ConditionList::iterator it=conditions.begin(); it!=conditions.end(); ++it){
		if((*it)->containsAtLeasetOneParameter()){
			return true;
		}
	}
	return false;
}

Dependency::ParameterParentMap BinaryLogicalCondition::getAllParameters(){
	Dependency::ParameterParentMap toReturn;
	Dependency::ParameterParentMap currentMap;
	for(ConditionList::iterator it = conditions.begin(); it != conditions.end(); ++it){
		currentMap = (*it)->getAllParameters();
		toReturn.insert(currentMap.begin(), currentMap.end());
	}
	return toReturn;
}

OrCondition::OrCondition(ConditionList& conditions):
	BinaryLogicalCondition(conditions){}

bool OrCondition::isConditionTrue(){
	ConditionList::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn || (*it)->isConditionTrue());
	}
	return toReturn;
}

AndCondition::AndCondition(ConditionList& conditions):
	BinaryLogicalCondition(conditions){}

bool AndCondition::isConditionTrue(){
	ConditionList::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn && (*it)->isConditionTrue());
	}
	return toReturn;
}

EqualsCondition::EqualsCondition(ConditionList& conditions):
	BinaryLogicalCondition(conditions){}

bool EqualsCondition::isConditionTrue(){
	ConditionList::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn == (*it)->isConditionTrue());
	}
	return toReturn;
}

NotCondition::NotCondition(Teuchos::RCP<Condition> condition):
	Condition(Condition::NotCon),
	condition(condition){}

bool NotCondition::isConditionTrue(){
	return (!condition->isConditionTrue());
}

bool NotCondition::containsAtLeasetOneParameter(){
	return condition->getType() == Condition::ParamCon;
}

Dependency::ParameterParentMap NotCondition::getAllParameters(){
	return condition->getAllParameters();
}

StringCondition::StringCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, std::string value, bool whenParamEqualsValue):
	ParameterCondition(parameterName, parentList, whenParamEqualsValue), values(ValueList(1,value))
{
	if(!parameter->isType<std::string>()){
		throw InvalidConditionException("The parameter of a String Condition "
		"must be of type string. \n"
		"Expected type: std::string\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

StringCondition::StringCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, ValueList values, bool whenParamEqualsValue):
	ParameterCondition(parameterName, parentList, whenParamEqualsValue), values(values)
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

BoolCondition::BoolCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue):
	ParameterCondition(parameterName, parentList, whenParamEqualsValue)
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

