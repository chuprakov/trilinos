
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


#ifndef CAZTECOO_STATUSTEST_CPP_HPP
#define CAZTECOO_STATUSTEST_CPP_HPP


#include "CTrilinos_config.h"

#ifdef HAVE_CTRILINOS_AZTECOO


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO_StatusTest.h"


namespace CAztecOO {


using Teuchos::RCP;


/*! get AztecOO_StatusTest from non-const table using CT_AztecOO_StatusTest_ID */
const RCP<AztecOO_StatusTest>
getStatusTest( CT_AztecOO_StatusTest_ID_t id );

/*! get AztecOO_StatusTest from non-const table using CTrilinos_Universal_ID_t */
const RCP<AztecOO_StatusTest>
getStatusTest( CTrilinos_Universal_ID_t id );

/*! get const AztecOO_StatusTest from either the const or non-const table
 * using CT_AztecOO_StatusTest_ID */
const RCP<const AztecOO_StatusTest>
getConstStatusTest( CT_AztecOO_StatusTest_ID_t id );

/*! get const AztecOO_StatusTest from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const AztecOO_StatusTest>
getConstStatusTest( CTrilinos_Universal_ID_t id );

/*! store AztecOO_StatusTest (owned) in non-const table */
CT_AztecOO_StatusTest_ID_t
storeNewStatusTest( AztecOO_StatusTest *pobj );

/*! store AztecOO_StatusTest in non-const table */
CT_AztecOO_StatusTest_ID_t
storeStatusTest( AztecOO_StatusTest *pobj );

/*! store const AztecOO_StatusTest in const table */
CT_AztecOO_StatusTest_ID_t
storeConstStatusTest( const AztecOO_StatusTest *pobj );

/* remove AztecOO_StatusTest from table using CT_AztecOO_StatusTest_ID */
void
removeStatusTest( CT_AztecOO_StatusTest_ID_t *id );

/* remove AztecOO_StatusTest from table using CTrilinos_Universal_ID_t */
void
removeStatusTest( CTrilinos_Universal_ID_t *aid );

/* purge AztecOO_StatusTest table */
void
purgeStatusTest(  );

/* store AztecOO_StatusTest in non-const table */
CTrilinos_Universal_ID_t
aliasStatusTest( const Teuchos::RCP< AztecOO_StatusTest > & robj );

/* store const AztecOO_StatusTest in const table */
CTrilinos_Universal_ID_t
aliasConstStatusTest( const Teuchos::RCP< const AztecOO_StatusTest > & robj );

} // namespace CAztecOO


#endif /* HAVE_CTRILINOS_AZTECOO */

#endif // CAZTECOO_STATUSTEST_CPP_HPP


