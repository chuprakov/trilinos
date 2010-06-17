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
#ifndef OPTIKA_STANDARDCONDITION_HPP_
#define OPTIKA_STANDARDCONDITION_HPP_
#include "Optika_Condition.hpp"
#include "Optika_InvalidConditionException.hpp"
#include "Teuchos_ParameterEntry.hpp"

namespace Optika{

class ParameterCondition : public Condition{
public:
	ParameterCondition(Teuchos::ParameterEntry* parameter, bool whenParamEqualsValue);

	virtual ~ParameterCondition(){}

	virtual bool isConditionTrue()=0;

protected:
	Teuchos::ParameterEntry* parameter;

	bool whenParamEqualsValue;
};

class BinaryLogicalCondition : public Condition{
public:
	BinaryLogicalCondition(std::vector<Teuchos::RCP<Condition> >& conditions);

	virtual ~BinaryLogicalCondition(){}
	
	virtual void addCondition(Teuchos::RCP<Condition> toAdd);

	virtual bool isConditionTrue()=0;

protected:
	std::vector<Teuchos::RCP<Condition> > conditions;
};

class OrCondition : public BinaryLogicalCondition{
public:
	OrCondition(std::vector<Teuchos::RCP<Condition> >& conditions);

	virtual ~OrCondition(){}

	bool isConditionTrue();
};

class AndCondition : public BinaryLogicalCondition{
public:
	AndCondition(std::vector<Teuchos::RCP<Condition> >& conditions);

	virtual ~AndCondition(){}

	bool isConditionTrue();
};

class EqualsCondition : public BinaryLogicalCondition{
public:
	EqualsCondition(std::vector<Teuchos::RCP<Condition> >& conditions);

	virtual ~EqualsCondition(){}

	bool isConditionTrue();
};

class NotCondition : public Condition{
public:
	NotCondition(Teuchos::RCP<Condition> condition);

	virtual ~NotCondition(){}

	bool isConditionTrue();

private:
	Teuchos::RCP<Condition> condition;
};

class StringCondition : public ParameterCondition{
public:
	/**
	 * Convience typedef
	 */
	typedef Teuchos::Array<std::string> ValueList; 

	StringCondition(Teuchos::ParameterEntry* parameter, std::string value, bool whenParamEqualsValue=true);

	StringCondition(Teuchos::ParameterEntry* parameter, ValueList values, bool whenParamEqualsValue=true);

	virtual ~StringCondition(){}

	bool isConditionTrue();

private:
	ValueList values;
};

template<class T>
class NumberCondition : public ParameterCondition{
public:
	NumberCondition(Teuchos::ParameterEntry* parameter, T (*func)(T)=0, bool whenParamEqualsValue=true):
		ParameterCondition(parameter, whenParamEqualsValue), func(func)
	{
		if(!parameter->isType<int>()
		&& !parameter->isType<short>()
		&& !parameter->isType<double>()
		&& !parameter->isType<float>()){
			throw InvalidConditionException("The parameter of a "
			"Number Condition must be of a supported number type!\n"
			"Expected type: " + typeid(T).name() + "\n" 
			"Actual type: " + parameter->getAny().typeName());
		}
	}
		

	NumberCondition(Teuchos::ParameterEntry* parameter, bool whenParamEqualsValue=true):
		ParameterCondition(parameter, whenParamEqualsValue), func(0){}

	virtual ~NumberCondition(){}

	bool isConditionTrue(){
		bool toReturn = (runFunction(Teuchos::getValue<T>(parameter)) > 0);
		return whenParamEqualsValue ? toReturn : !toReturn;
	}

private:
	T (*func)(T); 	

	T runFunction(T argument) const{
		if(func !=0)
			return (*func)(argument);
		else
			return argument;
	}	
};

class BoolCondition : public ParameterCondition{
public:
	BoolCondition(Teuchos::ParameterEntry* parameter, bool whenParamEqualsValue=true);

	virtual ~BoolCondition(){}

	bool isConditionTrue();
};

}
#endif //OPTIKA_STANDARDCONDITION_HPP_
