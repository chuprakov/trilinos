
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

#ifndef CEPETRA_CRSMATRIX_CPP_HPP
#define CEPETRA_CRSMATRIX_CPP_HPP


#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const RCP<Epetra_CrsMatrix>
getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id );

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_CrsMatrix>
getCrsMatrix( CTrilinos_Object_ID_t id );

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const RCP<const Epetra_CrsMatrix>
getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id );

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_CrsMatrix>
getConstCrsMatrix( CTrilinos_Object_ID_t id );

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
storeCrsMatrix( Epetra_CrsMatrix *pobj );

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
storeConstCrsMatrix( const Epetra_CrsMatrix *pobj );

/* dump contents of Epetra_CrsMatrix and const Epetra_CrsMatrix tables */
void
purgeCrsMatrixTables(  );


} // namespace CEpetra


#endif // CEPETRA_CRSMATRIX_CPP_HPP


