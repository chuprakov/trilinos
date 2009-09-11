
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

#ifndef CEPETRA_DISTRIBUTOR_CPP_HPP
#define CEPETRA_DISTRIBUTOR_CPP_HPP


#include "Epetra_Distributor.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_Distributor from non-const table using CT_Epetra_Distributor_ID */
const RCP<Epetra_Distributor>
getDistributor( CT_Epetra_Distributor_ID_t id );

/* get Epetra_Distributor from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_Distributor>
getDistributor( CTrilinos_Object_ID_t id );

/* get const Epetra_Distributor from either the const or non-const table
 * using CT_Epetra_Distributor_ID */
const RCP<const Epetra_Distributor>
getConstDistributor( CT_Epetra_Distributor_ID_t id );

/* get const Epetra_Distributor from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_Distributor>
getConstDistributor( CTrilinos_Object_ID_t id );

/* store Epetra_Distributor in non-const table */
CT_Epetra_Distributor_ID_t
storeDistributor( Epetra_Distributor *pobj );

/* store const Epetra_Distributor in const table */
CT_Epetra_Distributor_ID_t
storeConstDistributor( const Epetra_Distributor *pobj );

/* dump contents of Epetra_Distributor and const Epetra_Distributor tables */
void
purgeDistributorTables(  );


} // namespace CEpetra


#endif // CEPETRA_DISTRIBUTOR_CPP_HPP


