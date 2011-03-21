
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


#ifndef CEPETRA_SERIALDENSEVECTOR_CPP_HPP
#define CEPETRA_SERIALDENSEVECTOR_CPP_HPP


#include "CTrilinos_config.h"


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_SerialDenseVector from non-const table using CT_Epetra_SerialDenseVector_ID */
const RCP<Epetra_SerialDenseVector>
getSerialDenseVector( CT_Epetra_SerialDenseVector_ID_t id );

/*! get Epetra_SerialDenseVector from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_SerialDenseVector>
getSerialDenseVector( CTrilinos_Universal_ID_t id );

/*! get const Epetra_SerialDenseVector from either the const or non-const table
 * using CT_Epetra_SerialDenseVector_ID */
const RCP<const Epetra_SerialDenseVector>
getConstSerialDenseVector( CT_Epetra_SerialDenseVector_ID_t id );

/*! get const Epetra_SerialDenseVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_SerialDenseVector>
getConstSerialDenseVector( CTrilinos_Universal_ID_t id );

/*! store Epetra_SerialDenseVector (owned) in non-const table */
CT_Epetra_SerialDenseVector_ID_t
storeNewSerialDenseVector( Epetra_SerialDenseVector *pobj );

/*! store Epetra_SerialDenseVector in non-const table */
CT_Epetra_SerialDenseVector_ID_t
storeSerialDenseVector( Epetra_SerialDenseVector *pobj );

/*! store const Epetra_SerialDenseVector in const table */
CT_Epetra_SerialDenseVector_ID_t
storeConstSerialDenseVector( const Epetra_SerialDenseVector *pobj );

/* remove Epetra_SerialDenseVector from table using CT_Epetra_SerialDenseVector_ID */
void
removeSerialDenseVector( CT_Epetra_SerialDenseVector_ID_t *id );

/* remove Epetra_SerialDenseVector from table using CTrilinos_Universal_ID_t */
void
removeSerialDenseVector( CTrilinos_Universal_ID_t *aid );

/* purge Epetra_SerialDenseVector table */
void
purgeSerialDenseVector(  );

/* store Epetra_SerialDenseVector in non-const table */
CTrilinos_Universal_ID_t
aliasSerialDenseVector( const Teuchos::RCP< Epetra_SerialDenseVector > & robj );

/* store const Epetra_SerialDenseVector in const table */
CTrilinos_Universal_ID_t
aliasConstSerialDenseVector( const Teuchos::RCP< const Epetra_SerialDenseVector > & robj );

} // namespace CEpetra


#endif // CEPETRA_SERIALDENSEVECTOR_CPP_HPP


