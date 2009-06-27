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


Table<Epetra_SerialComm>& tableOfSerialComms()
{
    static Table<Epetra_SerialComm> loc_tableOfSerialComms(CT_Epetra_SerialComm_ID, "CT_Epetra_SerialComm_ID");
    return loc_tableOfSerialComms;
}


} // namespace


//
// Definitions from CEpetra_SerialComm.h
//


extern "C" {


CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfSerialComms(), id);
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  )
{
    return tableOfSerialComms().store(new Epetra_SerialComm());

}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID )
{
    const Teuchos::RCP<Epetra_SerialComm> 
        pComm = CEpetra::getSerialComm(CommID);

    return tableOfSerialComms().store(new Epetra_SerialComm(*pComm));

}

CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeComm(
        CEpetra::getSerialComm(selfID)->Clone());
}

void Epetra_SerialComm_Destroy ( 
  CT_Epetra_SerialComm_ID_t * selfID )
{
    tableOfSerialComms().remove(selfID);
}

void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID )
{
    CEpetra::getSerialComm(selfID)->Barrier();
}

int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root )
{
    return CEpetra::getSerialComm(selfID)->Broadcast(
        MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root )
{
    return CEpetra::getSerialComm(selfID)->Broadcast(
        MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root )
{
    return CEpetra::getSerialComm(selfID)->Broadcast(
        MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root )
{
    return CEpetra::getSerialComm(selfID)->Broadcast(
        MyVals, Count, Root);
}

int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getSerialComm(selfID)->GatherAll(
        MyVals, AllVals, Count);
}

int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getSerialComm(selfID)->GatherAll(
        MyVals, AllVals, Count);
}

int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getSerialComm(selfID)->GatherAll(
        MyVals, AllVals, Count);
}

int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getSerialComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count )
{
    return CEpetra::getSerialComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count )
{
    return CEpetra::getSerialComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getSerialComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count )
{
    return CEpetra::getSerialComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count )
{
    return CEpetra::getSerialComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getSerialComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count )
{
    return CEpetra::getSerialComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count )
{
    return CEpetra::getSerialComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count )
{
    return CEpetra::getSerialComm(selfID)->ScanSum(
        MyVals, ScanSums, Count);
}

int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getSerialComm(selfID)->ScanSum(
        MyVals, ScanSums, Count);
}

int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getSerialComm(selfID)->ScanSum(
        MyVals, ScanSums, Count);
}

int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::getSerialComm(selfID)->MyPID();
}

int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::getSerialComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeDistributor(
        CEpetra::getSerialComm(selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return CEpetra::storeDirectory(
        CEpetra::getSerialComm(selfID)->CreateDirectory(
            *pMap));
}

void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_SerialComm_ID_t CommID )
{
    Epetra_SerialComm& self = *( CEpetra::getSerialComm(selfID) );

    const Teuchos::RCP<Epetra_SerialComm> 
        pComm = CEpetra::getSerialComm(CommID);

    self = *pComm;
}


} // extern "C"


//
// Definitions from CEpetra_SerialComm_Cpp.hpp
//


const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return tableOfSerialComms().get(id);
}

CT_Epetra_SerialComm_ID_t
CEpetra::storeSerialComm( const Epetra_SerialComm *pobj )
{
    return tableOfSerialComms().storeCopy(pobj);
}

void
CEpetra::purgeSerialCommTable(  )
{
    tableOfSerialComms().purge();
}



