#include "CTrilinos_config.h"

#include "CTeuchos_any_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_ParameterEntry.h"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Teuchos::ParameterEntry */
Table<Teuchos::ParameterEntry>& tableOfParameterEntrys()
{
    static Table<Teuchos::ParameterEntry>
        loc_tableOfParameterEntrys(CT_Teuchos_ParameterEntry_ID, "CT_Teuchos_ParameterEntry_ID", false);
    return loc_tableOfParameterEntrys;
}

/* table to hold objects of type const Teuchos::ParameterEntry */
Table<const Teuchos::ParameterEntry>& tableOfConstParameterEntrys()
{
    static Table<const Teuchos::ParameterEntry>
        loc_tableOfConstParameterEntrys(CT_Teuchos_ParameterEntry_ID, "CT_Teuchos_ParameterEntry_ID", true);
    return loc_tableOfConstParameterEntrys;
}


} // namespace


//
// Definitions from CTeuchos_ParameterEntry.h
//


extern "C" {


CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterEntry_ID_t>(
        tableOfParameterEntrys().store(
        new Teuchos::ParameterEntry()));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( 
  CT_Teuchos_ParameterEntry_ID_t sourceID )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterEntry_ID_t>(
        tableOfParameterEntrys().store(new Teuchos::ParameterEntry(
        *CTeuchos::getConstParameterEntry(sourceID))));
}

void Teuchos_ParameterEntry_Destroy ( 
  CT_Teuchos_ParameterEntry_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Teuchos_ParameterEntry_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstParameterEntrys().remove(&aid);
    } else {
        tableOfParameterEntrys().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Teuchos_ParameterEntry_ID_t>(aid);
}

void Teuchos_ParameterEntry_Assign ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_ParameterEntry_ID_t sourceID )
{
    Teuchos::ParameterEntry& self = *( CTeuchos::getParameterEntry(selfID) );

    self = *CTeuchos::getConstParameterEntry(sourceID);
}

void Teuchos_ParameterEntry_setAnyValue ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_any_ID_t valueID, boolean isDefault )
{
    CTeuchos::getParameterEntry(selfID)->setAnyValue(
        *CTeuchos::getConstany(valueID), isDefault);
}

void Teuchos_ParameterEntry_setDocString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, char docString[] )
{
    CTeuchos::getParameterEntry(selfID)->setDocString(std::string(
        docString));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean isDefault, 
  char docString[] )
{
    return CTeuchos::storeParameterList(
        &( CTeuchos::getParameterEntry(selfID)->setList(isDefault, 
        std::string(docString)) ));
}

double Teuchos_ParameterEntry_getValue_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, double * ptr )
{
    return CTeuchos::getConstParameterEntry(selfID)->getValue<double>(
        ptr);
}

int Teuchos_ParameterEntry_getValue_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, int * ptr )
{
    return CTeuchos::getConstParameterEntry(selfID)->getValue<int>(
        ptr);
}

CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry )
{
    return CTeuchos::storeany(&( CTeuchos::getParameterEntry(
        selfID)->getAny(activeQry) ));
}

CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry )
{
    return CTeuchos::storeConstany(
        &( CTeuchos::getConstParameterEntry(selfID)->getAny(
        activeQry) ));
}

boolean Teuchos_ParameterEntry_isUsed ( 
  CT_Teuchos_ParameterEntry_ID_t selfID )
{
    return CTeuchos::getConstParameterEntry(selfID)->isUsed();
}

boolean Teuchos_ParameterEntry_isList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID )
{
    return CTeuchos::getConstParameterEntry(selfID)->isList();
}

boolean Teuchos_ParameterEntry_isType_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID )
{
    return CTeuchos::getConstParameterEntry(
        selfID)->isType<double>();
}

boolean Teuchos_ParameterEntry_isType_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID )
{
    return CTeuchos::getConstParameterEntry(selfID)->isType<int>();
}

boolean Teuchos_ParameterEntry_isDefault ( 
  CT_Teuchos_ParameterEntry_ID_t selfID )
{
    return CTeuchos::getConstParameterEntry(selfID)->isDefault();
}

const char * Teuchos_ParameterEntry_docString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID )
{
    return CTeuchos::getConstParameterEntry(
        selfID)->docString().c_str();
}


} // extern "C"


//
// Definitions from CTeuchos_ParameterEntry_Cpp.hpp
//


/* get Teuchos::ParameterEntry from non-const table using CT_Teuchos_ParameterEntry_ID */
const Teuchos::RCP<Teuchos::ParameterEntry>
CTeuchos::getParameterEntry( CT_Teuchos_ParameterEntry_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Teuchos_ParameterEntry_ID_t>(id);
    return tableOfParameterEntrys().get(aid);
}

/* get Teuchos::ParameterEntry from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Teuchos::ParameterEntry>
CTeuchos::getParameterEntry( CTrilinos_Object_ID_t id )
{
    return tableOfParameterEntrys().get(id);
}

/* get const Teuchos::ParameterEntry from either the const or non-const table
 * using CT_Teuchos_ParameterEntry_ID */
const Teuchos::RCP<const Teuchos::ParameterEntry>
CTeuchos::getConstParameterEntry( CT_Teuchos_ParameterEntry_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Teuchos_ParameterEntry_ID_t>(id);
    if (id.is_const) {
        return tableOfConstParameterEntrys().get(aid);
    } else {
        return tableOfParameterEntrys().get(aid);
    }
}

/* get const Teuchos::ParameterEntry from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Teuchos::ParameterEntry>
CTeuchos::getConstParameterEntry( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstParameterEntrys().get(id);
    } else {
        return tableOfParameterEntrys().get(id);
    }
}

/* store Teuchos::ParameterEntry in non-const table */
CT_Teuchos_ParameterEntry_ID_t
CTeuchos::storeParameterEntry( Teuchos::ParameterEntry *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterEntry_ID_t>(
            tableOfParameterEntrys().storeShared(pobj));
}

/* store const Teuchos::ParameterEntry in const table */
CT_Teuchos_ParameterEntry_ID_t
CTeuchos::storeConstParameterEntry( const Teuchos::ParameterEntry *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterEntry_ID_t>(
            tableOfConstParameterEntrys().storeShared(pobj));
}

/* dump contents of Teuchos::ParameterEntry and const Teuchos::ParameterEntry tables */
void
CTeuchos::purgeParameterEntryTables(  )
{
    tableOfParameterEntrys().purge();
    tableOfConstParameterEntrys().purge();
}



