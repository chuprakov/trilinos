
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


#ifdef HAVE_MPI


#include "CTrilinos_enums.h"
#include "CEpetra_MpiComm.h"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "mpi.h"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


//
// Definitions from CEpetra_MpiComm.h
//


extern "C" {


CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_MpiComm_Generalize ( 
  CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id);
}

CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm )
{
    return CEpetra::storeNewMpiComm(new Epetra_MpiComm(comm));
}

CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_MpiComm> Comm = CEpetra::getConstMpiComm(
        CommID);
    return CEpetra::storeNewMpiComm(new Epetra_MpiComm(*Comm));
}

void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID )
{
    CEpetra::removeMpiComm(selfID);
}

CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( 
  CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::storeComm(CEpetra::getConstMpiComm(selfID)->Clone());
}

void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID )
{
    CEpetra::getConstMpiComm(selfID)->Barrier();
}

int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root )
{
    return CEpetra::getConstMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root )
{
    return CEpetra::getConstMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root )
{
    return CEpetra::getConstMpiComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getConstMpiComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getConstMpiComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getConstMpiComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count )
{
    return CEpetra::getConstMpiComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count )
{
    return CEpetra::getConstMpiComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getConstMpiComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getConstMpiComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getConstMpiComm(selfID)->Comm();
}

int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getConstMpiComm(selfID)->MyPID();
}

int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getConstMpiComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( 
  CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getConstMpiComm(
        selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeDirectory(CEpetra::getConstMpiComm(
        selfID)->CreateDirectory(*Map));
}

int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getConstMpiComm(selfID)->GetMpiTag();
}

MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID )
{
    return CEpetra::getConstMpiComm(selfID)->GetMpiComm();
}

void Epetra_MpiComm_Assign ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID )
{
    Epetra_MpiComm& self = *( CEpetra::getMpiComm(selfID) );

    const Teuchos::RCP<const Epetra_MpiComm> Comm = CEpetra::getConstMpiComm(
        CommID);
    self = *Comm;
}


} // extern "C"


//
// Definitions from CEpetra_MpiComm_Cpp.hpp
//


/* get Epetra_MpiComm from non-const table using CT_Epetra_MpiComm_ID */
const Teuchos::RCP<Epetra_MpiComm>
CEpetra::getMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MpiComm>(
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id));
}

/* get Epetra_MpiComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MpiComm>
CEpetra::getMpiComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MpiComm>(id);
}

/* get const Epetra_MpiComm from either the const or non-const table
 * using CT_Epetra_MpiComm_ID */
const Teuchos::RCP<const Epetra_MpiComm>
CEpetra::getConstMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MpiComm>(
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id));
}

/* get const Epetra_MpiComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MpiComm>
CEpetra::getConstMpiComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MpiComm>(id);
}

/* store Epetra_MpiComm (owned) in non-const table */
CT_Epetra_MpiComm_ID_t
CEpetra::storeNewMpiComm( Epetra_MpiComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MpiComm>(pobj, true));
}

/* store Epetra_MpiComm in non-const table */
CT_Epetra_MpiComm_ID_t
CEpetra::storeMpiComm( Epetra_MpiComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MpiComm>(pobj, false));
}

/* store const Epetra_MpiComm in const table */
CT_Epetra_MpiComm_ID_t
CEpetra::storeConstMpiComm( const Epetra_MpiComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MpiComm>(pobj, false));
}

/* remove Epetra_MpiComm from table using CT_Epetra_MpiComm_ID */
void
CEpetra::removeMpiComm( CT_Epetra_MpiComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(aid);
}

/* purge Epetra_MpiComm table */
void
CEpetra::purgeMpiComm(  )
{
    CTrilinos::tableRepos().purge<Epetra_MpiComm>();
}



#endif /* HAVE_MPI */


