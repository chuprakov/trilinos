
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


#ifdef HAVE_CTRILINOS_AMESOS


#ifndef CAMESOS_CPP_HPP
#define CAMESOS_CPP_HPP


#include "Amesos.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CAmesos {


using Teuchos::RCP;


/* get Amesos from non-const table using CT_Amesos_ID */
const RCP<Amesos>
getAmesos( CT_Amesos_ID_t id );

/* get Amesos from non-const table using CTrilinos_Object_ID_t */
const RCP<Amesos>
getAmesos( CTrilinos_Object_ID_t id );

/* get const Amesos from either the const or non-const table
 * using CT_Amesos_ID */
const RCP<const Amesos>
getConstAmesos( CT_Amesos_ID_t id );

/* get const Amesos from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Amesos>
getConstAmesos( CTrilinos_Object_ID_t id );

/* store Amesos in non-const table */
CT_Amesos_ID_t
storeAmesos( Amesos *pobj );

/* store const Amesos in const table */
CT_Amesos_ID_t
storeConstAmesos( const Amesos *pobj );

/* dump contents of Amesos and const Amesos tables */
void
purgeAmesosTables(  );


} // namespace CAmesos


#endif // CAMESOS_CPP_HPP


#endif /* HAVE_CTRILINOS_AMESOS */


