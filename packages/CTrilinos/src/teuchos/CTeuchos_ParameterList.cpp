#include "CTrilinos_config.h"

#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CTeuchos_ParameterList.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Teuchos::ParameterList */
Table<Teuchos::ParameterList>& tableOfParameterLists()
{
    static Table<Teuchos::ParameterList>
        loc_tableOfParameterLists(CT_Teuchos_ParameterList_ID, "CT_Teuchos_ParameterList_ID", false);
    return loc_tableOfParameterLists;
}

/* table to hold objects of type const Teuchos::ParameterList */
Table<const Teuchos::ParameterList>& tableOfConstParameterLists()
{
    static Table<const Teuchos::ParameterList>
        loc_tableOfConstParameterLists(CT_Teuchos_ParameterList_ID, "CT_Teuchos_ParameterList_ID", true);
    return loc_tableOfConstParameterLists;
}


} // namespace


//
// Definitions from CTeuchos_ParameterList.h
//


extern "C" {


CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstParameterLists(), id);
    } else {
        newid = CTrilinos::cast(tableOfParameterLists(), id);
    }
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(newid);
}

CTrilinos_Object_ID_t Teuchos_ParameterList_Abstract ( 
  CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tableOfParameterLists().store(
        new Teuchos::ParameterList()));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tableOfParameterLists().store(new Teuchos::ParameterList(
        std::string(name))));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tableOfParameterLists().store(new Teuchos::ParameterList(
        *CTeuchos::getConstParameterList(sourceID))));
}

void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstParameterLists().remove(&aid);
    } else {
        tableOfParameterLists().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(aid);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setName(std::string(name)) ));
}

void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    Teuchos::ParameterList& self = *( CTeuchos::getParameterList(selfID) );

    self = *CTeuchos::getConstParameterList(sourceID);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setParameters(*CTeuchos::getConstParameterList(
        sourceID)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setParametersNotAlreadySet(
        *CTeuchos::getConstParameterList(sourceID)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->disableRecursiveValidation() ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->set(std::string(name), 
        *CTeuchos::getConstParameterList(valueID), std::string(
        docString)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setEntry(std::string(name), 
        *CTeuchos::getConstParameterEntry(entryID)) ));
}

double Teuchos_ParameterList_get_double_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value )
{
    return CTeuchos::getParameterList(selfID)->get<double>(
        std::string(name), def_value);
}

int Teuchos_ParameterList_get_int_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value )
{
    return CTeuchos::getParameterList(selfID)->get<int>(std::string(
        name), def_value);
}

const char * Teuchos_ParameterList_get_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] )
{
    return CTeuchos::getParameterList(selfID)->get(std::string(name), 
        def_value).c_str();
}

const char * Teuchos_ParameterList_get_const_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] )
{
    return CTeuchos::getParameterList(selfID)->get(std::string(name), 
        def_value).c_str();
}

double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->get<double>(
        std::string(name));
}

int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->get<int>(std::string(
        name));
}

double Teuchos_ParameterList_get_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->get<double>(
        std::string(name));
}

int Teuchos_ParameterList_get_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->get<int>(
        std::string(name));
}

double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->getPtr<double>(
        std::string(name));
}

int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->getPtr<int>(
        std::string(name));
}

const double * Teuchos_ParameterList_getPtr_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->getPtr<double>(
        std::string(name));
}

const int * Teuchos_ParameterList_getPtr_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->getPtr<int>(
        std::string(name));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeParameterEntry(
        &( CTeuchos::getParameterList(selfID)->getEntry(std::string(
        name)) ));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeConstParameterEntry(
        &( CTeuchos::getConstParameterList(selfID)->getEntry(
        std::string(name)) ));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeParameterEntry(CTeuchos::getParameterList(
        selfID)->getEntryPtr(std::string(name)));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeConstParameterEntry(
        CTeuchos::getConstParameterList(selfID)->getEntryPtr(
        std::string(name)));
}

boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists )
{
    return CTeuchos::getParameterList(selfID)->remove(std::string(
        name), throwIfNotExists);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->sublist(std::string(name), mustAlreadyExist, 
        std::string(docString)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeConstParameterList(
        &( CTeuchos::getConstParameterList(selfID)->sublist(
        std::string(name)) ));
}

const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTeuchos::getConstParameterList(selfID)->name().c_str();
}

boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->isParameter(
        std::string(name));
}

boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->isSublist(
        std::string(name));
}

boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->isType<double>(
        std::string(name));
}

boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->isType<int>(
        std::string(name));
}

boolean Teuchos_ParameterList_isType_double_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr )
{
    return CTeuchos::getConstParameterList(selfID)->isType<double>(
        std::string(name), ptr);
}

boolean Teuchos_ParameterList_isType_int_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr )
{
    return CTeuchos::getConstParameterList(selfID)->isType<int>(
        std::string(name), ptr);
}

const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTeuchos::getConstParameterList(
        selfID)->currentParametersString().c_str();
}

void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults )
{
    CTeuchos::getConstParameterList(selfID)->validateParameters(
        *CTeuchos::getConstParameterList(validParamListID), depth, 
        (Teuchos::EValidateUsed) validateUsed, 
        (Teuchos::EValidateDefaults) validateDefaults);
}

void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth )
{
    CTeuchos::getParameterList(
        selfID)->validateParametersAndSetDefaults(
        *CTeuchos::getConstParameterList(validParamListID), depth);
}


} // extern "C"


//
// Definitions from CTeuchos_ParameterList_Cpp.hpp
//


/* get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id);
    return tableOfParameterLists().get(aid);
}

/* get Teuchos::ParameterList from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CTrilinos_Object_ID_t id )
{
    return tableOfParameterLists().get(id);
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id);
    if (id.is_const) {
        return tableOfConstParameterLists().get(aid);
    } else {
        return tableOfParameterLists().get(aid);
    }
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstParameterLists().get(id);
    } else {
        return tableOfParameterLists().get(id);
    }
}

/* store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeParameterList( Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
            tableOfParameterLists().storeShared(pobj));
}

/* store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeConstParameterList( const Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
            tableOfConstParameterLists().storeShared(pobj));
}

/* dump contents of Teuchos::ParameterList and const Teuchos::ParameterList tables */
void
CTeuchos::purgeParameterListTables(  )
{
    tableOfParameterLists().purge();
    tableOfConstParameterLists().purge();
}



