#ifndef CEPETRA_SERIALCOMM_CPP_HPP
#define CEPETRA_SERIALCOMM_CPP_HPP

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
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialComm.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const RCP<Epetra_SerialComm>
getSerialComm( CT_Epetra_SerialComm_ID_t id );

/*! get Epetra_SerialComm from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_SerialComm>
getSerialComm( CTrilinos_Universal_ID_t id );

/*! get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const RCP<const Epetra_SerialComm>
getConstSerialComm( CT_Epetra_SerialComm_ID_t id );

/*! get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_SerialComm>
getConstSerialComm( CTrilinos_Universal_ID_t id );

/*! store Epetra_SerialComm (owned) in non-const table */
CT_Epetra_SerialComm_ID_t
storeNewSerialComm( Epetra_SerialComm *pobj );

/*! store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
storeSerialComm( Epetra_SerialComm *pobj );

/*! store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
storeConstSerialComm( const Epetra_SerialComm *pobj );

/* remove Epetra_SerialComm from table using CT_Epetra_SerialComm_ID */
void
removeSerialComm( CT_Epetra_SerialComm_ID_t *id );

/* remove Epetra_SerialComm from table using CTrilinos_Universal_ID_t */
void
removeSerialComm( CTrilinos_Universal_ID_t *aid );

/* purge Epetra_SerialComm table */
void
purgeSerialComm(  );

/* store Epetra_SerialComm in non-const table */
CTrilinos_Universal_ID_t
aliasSerialComm( const Teuchos::RCP< Epetra_SerialComm > & robj );

/* store const Epetra_SerialComm in const table */
CTrilinos_Universal_ID_t
aliasConstSerialComm( const Teuchos::RCP< const Epetra_SerialComm > & robj );

} // namespace CEpetra


#endif // CEPETRA_SERIALCOMM_CPP_HPP


