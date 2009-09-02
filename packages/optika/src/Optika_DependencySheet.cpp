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
#include "Optika_DependencySheet.hpp"
#include <QString>

namespace Optika{


DependencySheet::DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList):name("DEP_ANONYMOUS"), rootList(rootList){}

DependencySheet::DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList, const std::string &name):name(name), rootList(rootList){}

bool DependencySheet::addDependency(Teuchos::RCP<Optika::Dependency> dependency){
	validateExistanceInRoot(dependency);
	return dependencies[dependency->getDependee()].insert(dependency).second;
}

bool DependencySheet::removeDependency(Teuchos::RCP<Optika::Dependency> dependency){
	const Teuchos::ParameterEntry* dependee = dependency->getDependee();
	for(DepSet::iterator it = dependencies[dependee].begin(); it != dependencies[dependee].end(); it++){
		if((*it) == dependency){
			dependencies[dependee].erase(it);
			return true;
		}
	}
	return false;
}

bool DependencySheet::hasDependents(const Teuchos::ParameterEntry* dependee) const{
	return (dependencies.find(dependee) != dependencies.end() && dependencies.find(dependee)->second.size() > 0);
}

const DependencySheet::DepSet& DependencySheet::getDependenciesForParameter(const Teuchos::ParameterEntry* dependee) const{
	return dependencies.find(dependee)->second;
}

DependencySheet::DepMap::iterator DependencySheet::depBegin(){
	return dependencies.begin();
}

DependencySheet::DepMap::iterator DependencySheet::depEnd(){
	return dependencies.end();
}

DependencySheet::DepMap::const_iterator DependencySheet::depBegin() const{
	return dependencies.begin();
}

DependencySheet::DepMap::const_iterator DependencySheet::depEnd() const{
	return dependencies.end();
}

void DependencySheet::printDeps(){
	std::cout << "Dependency Sheet: " << name << "\n\n";
	for(DepMap::iterator it = depBegin(); it != depEnd(); it++){
		const Teuchos::ParameterEntry* dependee = it->first;
		for(DepSet::iterator it2 = dependencies.find(dependee)->second.begin(); it2 != dependencies.find(dependee)->second.end(); it2++){
			std::cout << "Dependee: " << (*it2)->getDependeeName() << "\n";
			std::cout << "Dependent: " << (*it2)->getDependentName() << "\n";
			std::cout << "Type: " << (*it2)->getType() << "\n\n";
		}
	}
}

void DependencySheet::validateExistanceInRoot(Teuchos::RCP<Optika::Dependency> dependency){
	if(!dependency->isDependeeParentInList(rootList)){
		throw InvalidDependencyException(
		"FAILED TO ADD DEPENDENCY!\n\n"
		"Sorry for the yelling there, but this is kind of a big deal. Dependencies are hard and complex so don't beat "
		"yourself up too much. Mistakes are easy to make when dealing with dependencies. "
		"And besides, I'm gonna do my best to help you out! I'm sure with the informationg below you'll be able to figure out what "
		"exactly went wrong. I've got confidence in you! :)\n\n"
		"Error:\n"
		"An attempt was made to add a dependency containing a the dependee parameter \"" + dependency->getDependeeName() + "\""
		" to the Dependency Sheet \"" + name + "\"."
		" The Dependecy Sheet's root list does not contain nor does it have"
		" child ParameterLists that contain the parameter.\n"
		"Dependency Sheet: " + name + "\n"
		"Dependency Type: " + QString::number(dependency->getType()).toStdString() + "\n"
		"Bad Parameter Name: " + dependency->getDependeeName());
	}
	if(!dependency->isDependentParentInList(rootList)){
		throw InvalidDependencyException(
		"FAILED TO ADD DEPENDENCY!\n\n"
		"Sorry for the yelling there, but this is kind of a big deal. Dependencies are hard and complex so don't beat "
		"yourself up too much. Mistakes are easy to make when dealing with dependencies. "
		"And besides, I'm gonna do my best to help you out! I'm sure with the informationg below you'll be able to figure out what "
		"exactly went wrong. I've got confidence in you! :)\n\n"
		"Error:\n"
		"An attempt was made to add a dependency containing a the dependent parameter \"" + dependency->getDependentName() + "\""
		" to the Dependency Sheet \"" + name + "\"."
		" The Dependency Sheet's list does not contain nor does it have"
		" child ParameterLists that contain the parameter.\n"
		"Dependency Sheet: " + name + "\n"
		"Dependency Type: " + QString::number(dependency->getType()).toStdString() + "\n"
		"Bad Parameter Name: " + dependency->getDependeeName());
	}
}



}

