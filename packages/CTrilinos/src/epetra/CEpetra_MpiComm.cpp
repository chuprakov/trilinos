#ifdef HAVE_MPI


#include "mpi.h"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_MpiComm.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_MpiComm>& tableOfMpiComms()
{
    static Table<Epetra_MpiComm> loc_tableOfMpiComms(CT_Epetra_MpiComm_ID, "CT_Epetra_MpiComm_ID");
    return loc_tableOfMpiComms;
}


} // namespace


//
// Definitions from CEpetra_MpiComm.h
//


extern "C" {


CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfMpiComms(), id);
}

CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm )
{
    return tableOfMpiComms().store(new Epetra_MpiComm(comm));

}

CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID )
{
    const Teuchos::RCP<Epetra_MpiComm> 
        pComm = CEpetra::getMpiComm(CommID);

    return tableOfMpiComms().store(new Epetra_MpiComm(*pComm));

}

void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID )
{
    tableOfMpiComms().remove(selfID);
}

void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID )
{
    CEpetra::getMpiComm(selfID)->Barrier();
}

int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root )
{
    return CEpetra::getMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root )
{
    return CEpetra::getMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root )
{
    return CEpetra::getMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root )
{
    return CEpetra::getMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->GatherAll(
        MyVals, AllVals, Count);
}

int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->GatherAll(
        MyVals, AllVals, Count);
}

int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->GatherAll(
        MyVals, AllVals, Count);
}

int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getMpiComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count )
{
    return CEpetra::getMpiComm(selfID)->SumAll(
        PartialSums, GlobalSums, Count);
}

int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getMpiComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count )
{
    return CEpetra::getMpiComm(selfID)->MaxAll(
        PartialMaxs, GlobalMaxs, Count);
}

int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getMpiComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count )
{
    return CEpetra::getMpiComm(selfID)->MinAll(
        PartialMins, GlobalMins, Count);
}

int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->ScanSum(
        MyVals, ScanSums, Count);
}

int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->ScanSum(
        MyVals, ScanSums, Count);
}

int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getMpiComm(selfID)->ScanSum(
        MyVals, ScanSums, Count);
}

MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getMpiComm(selfID)->Comm();
}

int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getMpiComm(selfID)->MyPID();
}

int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getMpiComm(selfID)->NumProc();
}

int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getMpiComm(selfID)->GetMpiTag();
}

MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getMpiComm(selfID)->GetMpiComm();
}


} // extern "C"


//
// Definitions from CEpetra_MpiComm_Cpp.hpp
//


const Teuchos::RCP<Epetra_MpiComm>
CEpetra::getMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return tableOfMpiComms().get(id);
}

void
CEpetra::purgeMpiCommTable(  )
{
    tableOfMpiComms().purge();
}



#endif /* HAVE_MPI */


