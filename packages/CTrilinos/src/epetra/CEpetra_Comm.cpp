
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

#include "CTrilinos_enums.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Comm_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


//
// Definitions from CEpetra_Comm.h
//


extern "C" {


CT_Epetra_Comm_ID_t Epetra_Comm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Comm_Generalize ( 
  CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id);
}

void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID )
{
    CEpetra::removeComm(selfID);
}

CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::storeComm(CEpetra::getConstComm(selfID)->Clone());
}

void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID )
{
    CEpetra::getConstComm(selfID)->Barrier();
}

int Epetra_Comm_Broadcast_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Char ( 
  CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_GatherAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_GatherAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_GatherAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_SumAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getConstComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_Comm_SumAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_Comm_SumAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialSums, long * GlobalSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_Comm_MaxAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getConstComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_Comm_MaxAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_Comm_MaxAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, long * GlobalMaxs, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_Comm_MinAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getConstComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_Comm_MinAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_Comm_MinAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMins, long * GlobalMins, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_Comm_ScanSum_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_ScanSum_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_ScanSum_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::getConstComm(selfID)->MyPID();
}

int Epetra_Comm_NumProc ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::getConstComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( 
  CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getConstComm(
        selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( 
  CT_Epetra_Comm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeDirectory(CEpetra::getConstComm(
        selfID)->CreateDirectory(*Map));
}


} // extern "C"


//
// Definitions from CEpetra_Comm_Cpp.hpp
//


/* get Epetra_Comm from non-const table using CT_Epetra_Comm_ID */
const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get Epetra_Comm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Comm>(id);
}

/* get const Epetra_Comm from either the const or non-const table
 * using CT_Epetra_Comm_ID */
const Teuchos::RCP<const Epetra_Comm>
CEpetra::getConstComm( CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get const Epetra_Comm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Comm>
CEpetra::getConstComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Comm>(id);
}

/* store Epetra_Comm (owned) in non-const table */
CT_Epetra_Comm_ID_t
CEpetra::storeNewComm( Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Comm>(pobj, true));
}

/* store Epetra_Comm in non-const table */
CT_Epetra_Comm_ID_t
CEpetra::storeComm( Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Comm>(pobj, false));
}

/* store const Epetra_Comm in const table */
CT_Epetra_Comm_ID_t
CEpetra::storeConstComm( const Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Comm>(pobj, false));
}

/* remove Epetra_Comm from table using CT_Epetra_Comm_ID */
void
CEpetra::removeComm( CT_Epetra_Comm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid);
}

/* purge Epetra_Comm table */
void
CEpetra::purgeComm(  )
{
    CTrilinos::tableRepos().purge<Epetra_Comm>();
}



