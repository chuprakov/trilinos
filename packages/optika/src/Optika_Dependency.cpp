// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Optika_Dependency.hpp"

namespace Optika{

Dependency::Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
		       std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, Type type)
{
	this->dependentParentList = dependentParentList;
	this->dependeeParentList = dependeeParentList;
	this->dependeeName = dependeeName;
	this->dependentName = dependentName;
	if(dependeeParentList->getEntryPtr(dependeeName) == NULL){
		throw InvalidDependencyException("The Dependee Parameter \"" + dependeeName + "\" does "
		"not exist in the given Dependee Parent List \"" + dependeeParentList->name() + "\"."
		"\n\nBummer! Maybe you just mispelled something? Why not go back and check to make sure "
		"you've got the names of the dependee and the depedent right? "
		"It might also be that you didn't specify the correct parent lists for the dependent and "
		"dependee. Either way, I'm sure it's just a simple mistake. "
		"You're a great programmer, I'm sure you'll figure it out! :)");
	}
	else{
		this->dependee = dependeeParentList->getEntryPtr(dependeeName);
	}
	if(dependentParentList->getEntryPtr(dependentName) == NULL){
		throw InvalidDependencyException("The Dependent Parameter \"" + dependentName + "\" does "
		"not exist in the given Dependent Parent List \"" + dependentParentList->name() + "\"."
		"\n\nBummer! Maybe you just mispelled something? Why not go back and check to make sure "
		"you've got the names of the dependee and the depedent right? "
		"It might also be that you just didn't specify the correct parent lists for the dependent and "
		"dependee. Either way, I'm sure it's just a simple mistake. "
		"You're a great programmer, I'm sure you'll figure it out! :)");
	}
	else{
		this->dependent = dependentParentList->getEntryPtr(dependentName);
	}
	this->type = type;
}

const Teuchos::ParameterEntry* Dependency::getDependee() const{
	return dependee;
}

const Teuchos::ParameterEntry* Dependency::getDependent() const{
	return dependent;
}

bool Dependency::isDependeeParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList){
	return doesListContainList(potentialParentList, dependeeParentList->name());
}

bool Dependency::isDependentParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList){
	return doesListContainList(potentialParentList, dependentParentList->name());
}

const std::string& Dependency::getDependeeName() const{
	return dependeeName;
}

const std::string& Dependency::getDependentName() const{
	return dependentName;
}

Dependency::Type Dependency::getType() const{
	return type;
}

bool Dependency::doesListContainList(Teuchos::RCP<Teuchos::ParameterList> parentList, std::string listName){
	if(parentList->name() == listName){
		return true;
	}
	else if(parentList->isSublist(listName)){
		return true;
	}
	else{
		for(Teuchos::ParameterList::ConstIterator it = parentList->begin(); it!=parentList->end(); it++){
			if(it->second.isList()){
				if(doesListContainList(sublist(parentList, it->first, true), listName)){
					return true;
				}
			}
		}
	}
	return false;
}



}

