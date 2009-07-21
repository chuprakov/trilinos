#include "CTrilinos_config.h"

#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_SerialComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_SerialComm */
Table<Epetra_SerialComm>& tableOfSerialComms()
{
    static Table<Epetra_SerialComm>
        loc_tableOfSerialComms(CT_Epetra_SerialComm_ID, "CT_Epetra_SerialComm_ID", false);
    return loc_tableOfSerialComms;
}

/* table to hold objects of type const Epetra_SerialComm */
Table<const Epetra_SerialComm>& tableOfConstSerialComms()
{
    static Table<const Epetra_SerialComm>
        loc_tableOfConstSerialComms(CT_Epetra_SerialComm_ID, "CT_Epetra_SerialComm_ID", true);
    return loc_tableOfConstSerialComms;
}


} // namespace


//
// Definitions from CEpetra_SerialComm.h
//


extern "C" {


CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstSerialComms(), id);
    } else {
        newid = CTrilinos::cast(tableOfSerialComms(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_SerialComm_Abstract ( 
  CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id);
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
            tableOfSerialComms().store(new Epetra_SerialComm()));
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
            tableOfSerialComms().store(new Epetra_SerialComm(
            *CEpetra::getSerialComm(CommID))));
}

CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeComm(
        CEpetra::getConstSerialComm(selfID)->Clone());
}

void Epetra_SerialComm_Destroy ( 
  CT_Epetra_SerialComm_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstSerialComms().remove(&aid);
    } else {
        tableOfSerialComms().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(aid);
}

void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID )
{
    CEpetra::getConstSerialComm(selfID)->Barrier();
}

int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(
            MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(
            MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(
            MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(
            MyVals, Count, Root);
}

int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(
            MyVals, AllVals, Count);
}

int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(
            MyVals, AllVals, Count);
}

int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(
            MyVals, AllVals, Count);
}

int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(
            PartialSums, GlobalSums, Count);
}

int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(
            PartialSums, GlobalSums, Count);
}

int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(
            PartialSums, GlobalSums, Count);
}

int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(
            PartialMaxs, GlobalMaxs, Count);
}

int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(
            PartialMaxs, GlobalMaxs, Count);
}

int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(
            PartialMaxs, GlobalMaxs, Count);
}

int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(
            PartialMins, GlobalMins, Count);
}

int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(
            PartialMins, GlobalMins, Count);
}

int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(
            PartialMins, GlobalMins, Count);
}

int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(
            MyVals, ScanSums, Count);
}

int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(
            MyVals, ScanSums, Count);
}

int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(
            MyVals, ScanSums, Count);
}

int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::getConstSerialComm(selfID)->MyPID();
}

int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::getConstSerialComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeDistributor(
        CEpetra::getConstSerialComm(selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    return CEpetra::storeDirectory(
        CEpetra::getConstSerialComm(selfID)->CreateDirectory(
                *CEpetra::getBlockMap(MapID)));
}

void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_SerialComm_ID_t CommID )
{
    Epetra_SerialComm& self = *( CEpetra::getSerialComm(selfID) );

    self = *CEpetra::getSerialComm(CommID);
}


} // extern "C"


//
// Definitions from CEpetra_SerialComm_Cpp.hpp
//


/* get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id);
    return tableOfSerialComms().get(aid);
}

/* get Epetra_SerialComm from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CTrilinos_Object_ID_t id )
{
    return tableOfSerialComms().get(id);
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id);
    if (id.is_const) {
        return tableOfConstSerialComms().get(aid);
    } else {
        return tableOfSerialComms().get(aid);
    }
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstSerialComms().get(id);
    } else {
        return tableOfSerialComms().get(id);
    }
}

/* store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
            tableOfSerialComms().storeCopy(pobj));
}

/* store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeConstSerialComm( const Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
            tableOfConstSerialComms().storeCopy(pobj));
}

/* dump contents of Epetra_SerialComm and const Epetra_SerialComm tables */
void
CEpetra::purgeSerialCommTables(  )
{
    tableOfSerialComms().purge();
    tableOfConstSerialComms().purge();
}



