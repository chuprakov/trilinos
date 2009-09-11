
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


#ifndef CEPETRA_MPICOMM_CPP_HPP
#define CEPETRA_MPICOMM_CPP_HPP


#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_MpiComm from non-const table using CT_Epetra_MpiComm_ID */
const RCP<Epetra_MpiComm>
getMpiComm( CT_Epetra_MpiComm_ID_t id );

/* get Epetra_MpiComm from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_MpiComm>
getMpiComm( CTrilinos_Object_ID_t id );

/* get const Epetra_MpiComm from either the const or non-const table
 * using CT_Epetra_MpiComm_ID */
const RCP<const Epetra_MpiComm>
getConstMpiComm( CT_Epetra_MpiComm_ID_t id );

/* get const Epetra_MpiComm from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_MpiComm>
getConstMpiComm( CTrilinos_Object_ID_t id );

/* store Epetra_MpiComm in non-const table */
CT_Epetra_MpiComm_ID_t
storeMpiComm( Epetra_MpiComm *pobj );

/* store const Epetra_MpiComm in const table */
CT_Epetra_MpiComm_ID_t
storeConstMpiComm( const Epetra_MpiComm *pobj );

/* dump contents of Epetra_MpiComm and const Epetra_MpiComm tables */
void
purgeMpiCommTables(  );


} // namespace CEpetra


#endif // CEPETRA_MPICOMM_CPP_HPP


#endif /* HAVE_MPI */


