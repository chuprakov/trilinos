
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


#ifdef HAVE_CTRILINOS_AZTECOO


#ifndef CAZTECOO_STATUSTESTMAXITERS_CPP_HPP
#define CAZTECOO_STATUSTESTMAXITERS_CPP_HPP


#include "AztecOO_StatusTestMaxIters.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CAztecOO {


using Teuchos::RCP;


/* get AztecOO_StatusTestMaxIters from non-const table using CT_AztecOO_StatusTestMaxIters_ID */
const RCP<AztecOO_StatusTestMaxIters>
getStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id );

/* get AztecOO_StatusTestMaxIters from non-const table using CTrilinos_Object_ID_t */
const RCP<AztecOO_StatusTestMaxIters>
getStatusTestMaxIters( CTrilinos_Object_ID_t id );

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CT_AztecOO_StatusTestMaxIters_ID */
const RCP<const AztecOO_StatusTestMaxIters>
getConstStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id );

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const AztecOO_StatusTestMaxIters>
getConstStatusTestMaxIters( CTrilinos_Object_ID_t id );

/* store AztecOO_StatusTestMaxIters in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
storeStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj );

/* store const AztecOO_StatusTestMaxIters in const table */
CT_AztecOO_StatusTestMaxIters_ID_t
storeConstStatusTestMaxIters( const AztecOO_StatusTestMaxIters *pobj );

/* dump contents of AztecOO_StatusTestMaxIters and const AztecOO_StatusTestMaxIters tables */
void
purgeStatusTestMaxItersTables(  );


} // namespace CAztecOO


#endif // CAZTECOO_STATUSTESTMAXITERS_CPP_HPP


#endif /* HAVE_CTRILINOS_AZTECOO */


