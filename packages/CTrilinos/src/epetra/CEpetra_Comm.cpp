#include "CTrilinos_config.h"

#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Comm.h"
#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Comm>& tableOfComms()
{
    static Table<Epetra_Comm> loc_tableOfComms(CT_Epetra_Comm_ID, "CT_Epetra_Comm_ID");
    return loc_tableOfComms;
}


} // namespace


//
// Definitions from CEpetra_Comm.h
//


extern "C" {


CT_Epetra_Comm_ID_t Epetra_Comm_Cast ( CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfComms(), id);
}

CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::storeComm(
        CEpetra::getComm(selfID)->Clone());
}

void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID )
{
    tableOfComms().remove(selfID);
}

void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID )
{
    CEpetra::getComm(selfID)->Barrier();
}

int Epetra_Comm_Broadcast_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, int Root )
{
    return CEpetra::getComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, int Root )
{
    return CEpetra::getComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, int Root )
{
    return CEpetra::getComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Char ( 
  CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, int Root )
{
    return CEpetra::getComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_GatherAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_GatherAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, int Count )
{
    return CEpetra::getComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_GatherAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_SumAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_Comm_SumAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count )
{
    return CEpetra::getComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_Comm_SumAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialSums, long * GlobalSums, 
  int Count )
{
    return CEpetra::getComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_Comm_MaxAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_Comm_MaxAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count )
{
    return CEpetra::getComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_Comm_MaxAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, long * GlobalMaxs, 
  int Count )
{
    return CEpetra::getComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_Comm_MinAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_Comm_MinAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count )
{
    return CEpetra::getComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_Comm_MinAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMins, long * GlobalMins, 
  int Count )
{
    return CEpetra::getComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_Comm_ScanSum_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count )
{
    return CEpetra::getComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_ScanSum_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_ScanSum_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::getComm(selfID)->MyPID();
}

int Epetra_Comm_NumProc ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::getComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( 
  CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::storeDistributor(
        CEpetra::getComm(selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( 
  CT_Epetra_Comm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return CEpetra::storeDirectory(
        CEpetra::getComm(selfID)->CreateDirectory(
            *pMap));
}


} // extern "C"


//
// Definitions from CEpetra_Comm_Cpp.hpp
//


const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CT_Epetra_Comm_ID_t id )
{
    return tableOfComms().get(id);
}

CT_Epetra_Comm_ID_t
CEpetra::storeComm( const Epetra_Comm *pobj )
{
    return tableOfComms().storeCopy(pobj);
}

void
CEpetra::purgeCommTable(  )
{
    tableOfComms().purge();
}



