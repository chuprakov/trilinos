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



#include "Teuchos_StandardConditions.hpp"

namespace Teuchos{

ParameterCondition::ParameterCondition(
  RCP<ParameterEntry> parameter, 
  bool whenParamEqualsValue):
  parameterEntry_(parameter),
  whenParamEqualsValue_(whenParamEqualsValue)
{
  TEST_FOR_EXCEPTION(is_null(parameter),
    InvalidConditionException,
    "Parameter conditions can't be given a null parameter" <<
    std::endl << std::endl);
}

Dependency::ConstParameterEntryList 
ParameterCondition::getAllParameters() const
{
  Dependency::ConstParameterEntryList toReturn;
  toReturn.insert(getParameter());
  return toReturn;
}

BinaryLogicalCondition::BinaryLogicalCondition(ConstConditionList& conditions):
  conditions_(conditions)
{
  if(conditions_.size() ==0){
    throw InvalidConditionException("Sorry bud, but you gotta at least give "
    "me one condition "
    "when you're constructing a BinaryLogicalCondition. Looks like you didn't. "
    "I'm just gonna "
    "chalk it up a silly little mistake though. Take a look over your "
    "conditions again and make sure "
    "you don't ever give any of your BinaryLogicConditions and empty "
    "condition list.\n\n"
    "Error: Empty condition list given to a BinaryLogicalCondition "
    "constructor.");
  }
}


void BinaryLogicalCondition::addCondition(RCP<const Condition> toAdd){
  conditions_.append(toAdd);
}

bool BinaryLogicalCondition::isConditionTrue() const{
  ConstConditionList::const_iterator it = conditions_.begin();
  bool toReturn = (*it)->isConditionTrue();
  ++it;
  for(;it != conditions_.end(); ++it){
    toReturn = applyOperator(toReturn,(*it)->isConditionTrue());
  }
  return toReturn;
}

bool BinaryLogicalCondition::containsAtLeasteOneParameter() const{
  for(
    ConstConditionList::const_iterator it=conditions_.begin();
    it!=conditions_.end();
    ++it)
  {
    if((*it)->containsAtLeasteOneParameter()){
      return true;
    }
  }
  return false;
}

Dependency::ConstParameterEntryList 
BinaryLogicalCondition::getAllParameters() const{
  Dependency::ConstParameterEntryList toReturn;
  Dependency::ConstParameterEntryList currentList;
  for(
    ConstConditionList::const_iterator it = conditions_.begin();
    it != conditions_.end();
    ++it)
  {
    currentList = (*it)->getAllParameters();
    toReturn.insert(currentList.begin(), currentList.end());
  }
  return toReturn;
}

OrCondition::OrCondition(ConstConditionList& conditions):
  BinaryLogicalCondition(conditions){}

bool OrCondition::applyOperator(bool op1, bool op2) const{
  return op1 || op2;
}

AndCondition::AndCondition(ConstConditionList& conditions):
  BinaryLogicalCondition(conditions){}

bool AndCondition::applyOperator(bool op1, bool op2) const{
  return op1 && op2;
}

EqualsCondition::EqualsCondition(ConstConditionList& conditions):
  BinaryLogicalCondition(conditions){}

bool EqualsCondition::applyOperator(bool op1, bool op2) const{
  return op1 == op2;
}

NotCondition::NotCondition(RCP<const Condition> childCondition):
  childCondition_(childCondition)
{
  if(childCondition_.is_null()){
    throw InvalidConditionException("OOOOOOOOPppppps! Looks like you tried "
    "to give me "
    "a null pointer when you were making a not conditon. "
    "That's a no no. Go back and "
    "checkout your not conditions and make sure you didn't "
    "give any of them a null pointer "
    "as an argument to the constructor.\n\n"
    "Error: Null pointer given to NotCondition constructor.");
  }
}

bool NotCondition::isConditionTrue() const{
  return (!childCondition_->isConditionTrue());
}

bool NotCondition::containsAtLeasteOneParameter() const{
  return childCondition_->containsAtLeasteOneParameter();
}

Dependency::ConstParameterEntryList NotCondition::getAllParameters() const{
  return childCondition_->getAllParameters();
}

StringCondition::StringCondition(
  RCP<ParameterEntry> parameter,
  std::string value, 
  bool whenParamEqualsValue):
  ParameterCondition(parameter, whenParamEqualsValue), 
  values_(ValueList(1,value))
{
  if(!getParameter()->isType<std::string>()){
    throw InvalidConditionException("The parameter of a String Condition "
    "must be of type string. \n"
    "Expected type: std::string\n"
    "Actual type: " + getParameter()->getAny().typeName() + "\n");
  }
}

StringCondition::StringCondition(
  RCP<ParameterEntry> parameter,
  ValueList values,
  bool whenParamEqualsValue):
  ParameterCondition(parameter, whenParamEqualsValue), 
  values_(values)
{
  if(!getParameter()->isType<std::string>()){
    throw InvalidConditionException("The parameter of a String Condition "
    "must be of type string. \n"
    "Expected type: std::string\n"
    "Actual type: " + getParameter()->getAny().typeName() + "\n");
  }
}

bool StringCondition::evaluateParameter() const{
  return find(
    values_.begin(), values_.end(), 
    getValue<std::string>(*getParameter())) != values_.end();
}

BoolCondition::BoolCondition(
  RCP<ParameterEntry> parameter,
  bool whenParamEqualsValue):
  ParameterCondition(parameter, whenParamEqualsValue)
{
  if(!getParameter()->isType<bool>()){
    throw InvalidConditionException("The parameter of a Bool Condition "
    "must be of type Bool. \n"
    "Expected type: Bool\n"
    "Actual type: " + getParameter()->getAny().typeName() + "\n");
  }
}

bool BoolCondition::evaluateParameter() const{
  return getValue<bool>(*getParameter());
}

} //namespace Teuchos

