
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
#include "CEpetra_Vector.h"
#include "CEpetra_Vector_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


//
// Definitions from CEpetra_Vector.h
//


extern "C" {


CT_Epetra_Vector_ID_t Epetra_Vector_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Vector_Generalize ( 
  CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id);
}

CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeNewVector(new Epetra_Vector(*Map, ((zeroOut) != 
        FALSE ? true : false)));
}

CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID )
{
    const Teuchos::RCP<const Epetra_Vector> Source = CEpetra::getConstVector(
        SourceID);
    return CEpetra::storeNewVector(new Epetra_Vector(*Source));
}

CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * V )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeNewVector(new Epetra_Vector((Epetra_DataAccess) CV, 
        *Map, V));
}

CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CEpetra::getConstMultiVector(SourceID);
    return CEpetra::storeNewVector(new Epetra_Vector((Epetra_DataAccess) CV, 
        *Source, Index));
}

void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID )
{
    CEpetra::removeVector(selfID);
}

int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceGlobalValues(NumEntries, Values, 
        Indices);
}

int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceMyValues(NumEntries, Values, 
        Indices);
}

int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoGlobalValues(NumEntries, Values, 
        Indices);
}

int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoMyValues(NumEntries, Values, 
        Indices);
}

int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceGlobalValues(NumEntries, 
        BlockOffset, Values, Indices);
}

int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceMyValues(NumEntries, BlockOffset, 
        Values, Indices);
}

int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoGlobalValues(NumEntries, 
        BlockOffset, Values, Indices);
}

int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoMyValues(NumEntries, BlockOffset, 
        Values, Indices);
}

int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V )
{
    return CEpetra::getConstVector(selfID)->ExtractCopy(V);
}

int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V )
{
    return CEpetra::getConstVector(selfID)->ExtractView(V);
}

double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index )
{
    const Epetra_Vector& self = *( CEpetra::getConstVector(selfID) );

    return self[index];
}


} // extern "C"


//
// Definitions from CEpetra_Vector_Cpp.hpp
//


/* get Epetra_Vector from non-const table using CT_Epetra_Vector_ID */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Vector>(
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id));
}

/* get Epetra_Vector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Vector>(id);
}

/* get const Epetra_Vector from either the const or non-const table
 * using CT_Epetra_Vector_ID */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Vector>(
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id));
}

/* get const Epetra_Vector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Vector>(id);
}

/* store Epetra_Vector (owned) in non-const table */
CT_Epetra_Vector_ID_t
CEpetra::storeNewVector( Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Vector>(pobj, true));
}

/* store Epetra_Vector in non-const table */
CT_Epetra_Vector_ID_t
CEpetra::storeVector( Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Vector>(pobj, false));
}

/* store const Epetra_Vector in const table */
CT_Epetra_Vector_ID_t
CEpetra::storeConstVector( const Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Vector>(pobj, false));
}

/* remove Epetra_Vector from table using CT_Epetra_Vector_ID */
void
CEpetra::removeVector( CT_Epetra_Vector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(aid);
}

/* purge Epetra_Vector table */
void
CEpetra::purgeVector(  )
{
    CTrilinos::tableRepos().purge<Epetra_Vector>();
}



