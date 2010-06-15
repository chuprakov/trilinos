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
#include "Optika_Dependency.hpp"

namespace Optika{

Dependency::Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap& dependents, Type type)/*:
	type(type),
	dependeeParentList(dependeeParentList),
	dependeeName(dependeeName)*/
{
	this->type = type;
	this->dependeeParentList = dependeeParentList;
	this->dependeeName = dependeeName;
	intitializeDependeesAndDependents(dependents);
}
	
Dependency::Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
 	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, Type type)/*:
	type(type),
	dependeeParentList(dependeeParentList),
	dependeeName(dependeeName)*/
{
	this->type = type;
	this->dependeeParentList = dependeeParentList;
	this->dependeeName = dependeeName;
	ParameterParentMap dependents;
	dependents.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(dependentName, dependentParentList));
	intitializeDependeesAndDependents(dependents);
}

void Dependency::intitializeDependeesAndDependents(ParameterParentMap& dependents){
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
	ParameterParentMap::iterator it;
	for(it = dependents.begin(); it != dependents.end(); ++it){
		if((*it).second->getEntryPtr((*it).first) == NULL){
			throw InvalidDependencyException("The Dependent Parameter \"" + (*it).first + "\" does "
			"not exist in the given Dependent Parent List \"" + (*it).second->name() + "\"."
			"\n\nBummer! Maybe you just mispelled something? Why not go back and check to make sure "
			"you've got the names of the dependee and the depedent right? "
			"It might also be that you just didn't specify the correct parent lists for the dependent and "
			"dependee. Either way, I'm sure it's just a simple mistake. "
			"You're a great programmer, I'm sure you'll figure it out! :)");
		}
		else{
			this->dependents.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(it->first, it->second));
			dependentNames.insert(it->first);
		}
	}
}


const Teuchos::ParameterEntry* Dependency::getDependee() const{
	return dependee;
}

Dependency::ParameterParentMap Dependency::getDependents() const{
	return dependents;
}

bool Dependency::isDependeeParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList){
	return doesListContainList(potentialParentList, dependeeParentList);
}

/*bool Dependency::areDependentsParentsInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList){
	ParameterParentMap::iterator it;
	for(it = dependents.begin(); it != dependents.end(); ++it){
		if(!doesListContainList(potentialParentList, it->second)){
			return false;
		}
	}
	return true;
}*/

const std::string& Dependency::getDependeeName() const{
	return dependeeName;
}

std::set<std::string> Dependency::getDependentNames() const{
	return dependentNames;
}

std::string Dependency::getDependentNamesString() const{
	std::string names = "";
	for(std::set<std::string>::const_iterator it=dependentNames.begin(); it != dependentNames.end(); ++it){
		names += *it + "" ;
	}
	return names;
}

Dependency::Type Dependency::getType() const{
	return type;
}

bool Dependency::doesListContainList(Teuchos::RCP<Teuchos::ParameterList> parentList, Teuchos::RCP<Teuchos::ParameterList> listToFind){
	if(parentList.get() == listToFind.get()){
		return true;
	}
	else{
		for(Teuchos::ParameterList::ConstIterator it = parentList->begin(); it!=parentList->end(); it++){
			if(it->second.isList()){
				if(doesListContainList(sublist(parentList, it->first,true), listToFind)){
					return true;
				}
			}
		}
	}
	return false;
}


}

