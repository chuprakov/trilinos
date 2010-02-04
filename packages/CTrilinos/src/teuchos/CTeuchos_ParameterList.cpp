
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"

#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CTeuchos_ParameterList.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CTeuchos_ParameterList.h
//


extern "C" {


CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(id);
}

CTrilinos_Universal_ID_t Teuchos_ParameterList_Generalize ( 
  CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(new Teuchos::ParameterList());
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(new Teuchos::ParameterList(std::string(
        name)));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(sourceID);
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(new Teuchos::ParameterList(*source));
}

void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->setName(std::string(name)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(sourceID);
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->setParameters(*source) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(sourceID);
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->setParametersNotAlreadySet(
        *source) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        selfID)->disableRecursiveValidation() ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  double value, char const docString[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->set<double>(std::string(name), 
        value, std::string(docString),  Teuchos::null ) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  int value, char const docString[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->set<int>(std::string(name), 
        value, std::string(docString),  Teuchos::null ) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_str ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  char value[], char const docString[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->set(std::string(name), value, 
        std::string(docString),  Teuchos::null ) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] )
{
    const Teuchos::RCP<const Teuchos::ParameterList> value = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(valueID);
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->set(std::string(name), *value, 
        std::string(docString)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID )
{
    const Teuchos::RCP<const Teuchos::ParameterEntry> entry = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterEntry, 
        CT_Teuchos_ParameterEntry_ID_t>(entryID);
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->setEntry(std::string(name), 
        *entry) ));
}

double Teuchos_ParameterList_get_double_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get<double>(std::string(name), 
        def_value);
}

int Teuchos_ParameterList_get_int_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get<int>(std::string(name), 
        def_value);
}

const char * Teuchos_ParameterList_get_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get(std::string(name), 
        def_value).c_str();
}

const char * Teuchos_ParameterList_get_const_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get(std::string(name), 
        def_value).c_str();
}

double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get<double>(std::string(
        name));
}

int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get<int>(std::string(name));
}

double Teuchos_ParameterList_get_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get<double>(std::string(
        name));
}

int Teuchos_ParameterList_get_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->get<int>(std::string(name));
}

double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getPtr<double>(std::string(
        name));
}

int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getPtr<int>(std::string(
        name));
}

const double * Teuchos_ParameterList_getPtr_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getPtr<double>(std::string(
        name));
}

const int * Teuchos_ParameterList_getPtr_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getPtr<int>(std::string(
        name));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterEntry, 
        CT_Teuchos_ParameterEntry_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getEntry(std::string(
        name)) ));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterEntry, 
        CT_Teuchos_ParameterEntry_ID_t>(
        &( CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getEntry(std::string(
        name)) ));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterEntry, 
        CT_Teuchos_ParameterEntry_ID_t>(
        CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getEntryPtr(std::string(
        name)));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterEntry, 
        CT_Teuchos_ParameterEntry_ID_t>(
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->getEntryPtr(std::string(
        name)));
}

boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists )
{
    return ((CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->remove(std::string(name), ((
        throwIfNotExists) != FALSE ? true : false))) ? TRUE : FALSE);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->sublist(std::string(name), ((
        mustAlreadyExist) != FALSE ? true : false), std::string(
        docString)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        &( CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->sublist(std::string(name)) ));
}

const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->name().c_str();
}

boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->isParameter(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->isSublist(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->isType<double>(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->isType<int>(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_double_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr )
{
    return ((CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->isType<double>(std::string(
        name), ptr)) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_int_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr )
{
    return ((CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->isType<int>(std::string(name), 
        ptr)) ? TRUE : FALSE);
}

const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        selfID)->currentParametersString().c_str();
}

void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults )
{
    const Teuchos::RCP<const Teuchos::ParameterList> validParamList = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(validParamListID);
    CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID)->validateParameters(
        *validParamList, depth, (Teuchos::EValidateUsed) validateUsed, 
        (Teuchos::EValidateDefaults) validateDefaults);
}

void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth )
{
    const Teuchos::RCP<const Teuchos::ParameterList> validParamList = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(validParamListID);
    CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(
        selfID)->validateParametersAndSetDefaults(*validParamList, depth);
}

void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    Teuchos::ParameterList& self = *( 
        CTrilinos::tableRepos().get<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(selfID) );

    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTrilinos::tableRepos().getConst<Teuchos::ParameterList, 
        CT_Teuchos_ParameterList_ID_t>(sourceID);
    self = *source;
}


} // extern "C"


//
// Definitions from CTeuchos_ParameterList_Cpp.hpp
//


/* get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, CT_Teuchos_ParameterList_ID_t>(id);
}

/* get Teuchos::ParameterList from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList, CTrilinos_Universal_ID_t>(id);
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, CT_Teuchos_ParameterList_ID_t>(id);
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList, CTrilinos_Universal_ID_t>(id);
}

/* store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeParameterList( Teuchos::ParameterList *pobj )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, CT_Teuchos_ParameterList_ID_t>(pobj, false);
}

/* store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeConstParameterList( const Teuchos::ParameterList *pobj )
{
    return CTrilinos::tableRepos().store<Teuchos::ParameterList, CT_Teuchos_ParameterList_ID_t>(pobj, false);
}



