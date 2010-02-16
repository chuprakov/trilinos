
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
#include "CEpetra_SerialComm.h"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


//
// Definitions from CEpetra_SerialComm.h
//


extern "C" {


CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_SerialComm_Generalize ( 
  CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id);
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  )
{
    return CEpetra::storeNewSerialComm(new Epetra_SerialComm());
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_SerialComm> Comm = 
        CEpetra::getConstSerialComm(CommID);
    return CEpetra::storeNewSerialComm(new Epetra_SerialComm(*Comm));
}

void Epetra_SerialComm_Destroy ( CT_Epetra_SerialComm_ID_t * selfID )
{
    CEpetra::removeSerialComm(selfID);
}

CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeComm(CEpetra::getConstSerialComm(selfID)->Clone());
}

void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID )
{
    CEpetra::getConstSerialComm(selfID)->Barrier();
}

int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * AllVals, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(MyVals, AllVals, 
        Count);
}

int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(MyVals, AllVals, 
        Count);
}

int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(MyVals, AllVals, 
        Count);
}

int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(MyVals, ScanSums, 
        Count);
}

int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(MyVals, ScanSums, 
        Count);
}

int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(MyVals, ScanSums, 
        Count);
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
    return CEpetra::storeDistributor(CEpetra::getConstSerialComm(
        selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeDirectory(CEpetra::getConstSerialComm(
        selfID)->CreateDirectory(*Map));
}

void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, 
  CT_Epetra_SerialComm_ID_t CommID )
{
    Epetra_SerialComm& self = *( CEpetra::getSerialComm(selfID) );

    const Teuchos::RCP<const Epetra_SerialComm> Comm = 
        CEpetra::getConstSerialComm(CommID);
    self = *Comm;
}


} // extern "C"


//
// Definitions from CEpetra_SerialComm_Cpp.hpp
//


/* get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get Epetra_SerialComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialComm>(id);
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SerialComm>(id);
}

/* store Epetra_SerialComm (owned) in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeNewSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialComm>(pobj, true));
}

/* store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialComm>(pobj, false));
}

/* store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeConstSerialComm( const Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialComm>(pobj, false));
}

/* remove Epetra_SerialComm from table using CT_Epetra_SerialComm_ID */
void
CEpetra::removeSerialComm( CT_Epetra_SerialComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(aid);
}

/* purge Epetra_SerialComm table */
void
CEpetra::purgeSerialComm(  )
{
    CTrilinos::tableRepos().purge<Epetra_SerialComm>();
}



