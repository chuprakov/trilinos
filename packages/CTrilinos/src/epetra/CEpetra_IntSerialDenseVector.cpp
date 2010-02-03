
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

#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_IntSerialDenseVector.h
//


extern "C" {


CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_Empty ( 
   )
{
    return CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector, 
        CT_Epetra_IntSerialDenseVector_ID_t>(
        new Epetra_IntSerialDenseVector());
}

CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create ( 
  int Length_in )
{
    return CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector, 
        CT_Epetra_IntSerialDenseVector_ID_t>(
        new Epetra_IntSerialDenseVector(Length_in));
}

CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV_in, int * Values_in, int Length_in )
{
    return CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector, 
        CT_Epetra_IntSerialDenseVector_ID_t>(
        new Epetra_IntSerialDenseVector((Epetra_DataAccess) CV_in, 
        Values_in, Length_in));
}

CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Duplicate ( 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID )
{
    return CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector, 
        CT_Epetra_IntSerialDenseVector_ID_t>(
        new Epetra_IntSerialDenseVector(
        *CEpetra::getConstIntSerialDenseVector(SourceID)));
}

void Epetra_IntSerialDenseVector_Destroy ( 
  CT_Epetra_IntSerialDenseVector_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_IntSerialDenseVector_Size ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Length_in )
{
    return CEpetra::getIntSerialDenseVector(selfID)->Size(Length_in);
}

int Epetra_IntSerialDenseVector_Resize ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Length_in )
{
    return CEpetra::getIntSerialDenseVector(selfID)->Resize(
        Length_in);
}

int Epetra_IntSerialDenseVector_Random ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID )
{
    return CEpetra::getIntSerialDenseVector(selfID)->Random();
}

int Epetra_IntSerialDenseVector_Length ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID )
{
    return CEpetra::getConstIntSerialDenseVector(selfID)->Length();
}

int * Epetra_IntSerialDenseVector_Values ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID )
{
    return CEpetra::getIntSerialDenseVector(selfID)->Values();
}

const int * Epetra_IntSerialDenseVector_Values_Const ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID )
{
    return CEpetra::getConstIntSerialDenseVector(selfID)->Values();
}

CT_Epetra_DataAccess_E_t Epetra_IntSerialDenseVector_CV ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID )
{
    return (CT_Epetra_DataAccess_E_t)(
         CEpetra::getConstIntSerialDenseVector(selfID)->CV() );
}

int Epetra_IntSerialDenseVector_MakeViewOf ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID )
{
    return CEpetra::getIntSerialDenseVector(selfID)->MakeViewOf(
        *CEpetra::getConstIntSerialDenseVector(SourceID));
}

void Epetra_IntSerialDenseVector_setElement ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index, 
  int * value )
{
    Epetra_IntSerialDenseVector& self = *( CEpetra::getIntSerialDenseVector(selfID) );

    self(Index) = *value;
}

int Epetra_IntSerialDenseVector_getElement ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index )
{
    const Epetra_IntSerialDenseVector& self = *( CEpetra::getConstIntSerialDenseVector(selfID) );

    return self(Index);
}

void Epetra_IntSerialDenseVector_setElement_Bracket ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index, 
  int * value )
{
    Epetra_IntSerialDenseVector& self = *( CEpetra::getIntSerialDenseVector(selfID) );

    self[Index] = *value;
}

int Epetra_IntSerialDenseVector_getElement_Bracket ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index )
{
    const Epetra_IntSerialDenseVector& self = *( CEpetra::getConstIntSerialDenseVector(selfID) );

    return self[Index];
}

void Epetra_IntSerialDenseVector_Assign ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID )
{
    Epetra_IntSerialDenseVector& self = *( CEpetra::getIntSerialDenseVector(selfID) );

    self = *CEpetra::getConstIntSerialDenseVector(SourceID);
}


} // extern "C"


//
// Definitions from CEpetra_IntSerialDenseVector_Cpp.hpp
//


/* get Epetra_IntSerialDenseVector from non-const table using CT_Epetra_IntSerialDenseVector_ID */
const Teuchos::RCP<Epetra_IntSerialDenseVector>
CEpetra::getIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_IntSerialDenseVector, CT_Epetra_IntSerialDenseVector_ID_t>(id);
}

/* get const Epetra_IntSerialDenseVector from either the const or non-const table
 * using CT_Epetra_IntSerialDenseVector_ID */
const Teuchos::RCP<const Epetra_IntSerialDenseVector>
CEpetra::getConstIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_IntSerialDenseVector, CT_Epetra_IntSerialDenseVector_ID_t>(id);
}

/* store Epetra_IntSerialDenseVector in non-const table */
CT_Epetra_IntSerialDenseVector_ID_t
CEpetra::storeIntSerialDenseVector( Epetra_IntSerialDenseVector *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector, CT_Epetra_IntSerialDenseVector_ID_t>(pobj, false);
}

/* store const Epetra_IntSerialDenseVector in const table */
CT_Epetra_IntSerialDenseVector_ID_t
CEpetra::storeConstIntSerialDenseVector( const Epetra_IntSerialDenseVector *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector, CT_Epetra_IntSerialDenseVector_ID_t>(pobj, false);
}



