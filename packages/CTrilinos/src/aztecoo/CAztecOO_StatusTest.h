
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



/*! @file CAztecOO_StatusTest.h
 * @brief Wrappers for AztecOO_StatusTest */

/* True C header file! */


#ifndef CAZTECOO_STATUSTEST_H
#define CAZTECOO_STATUSTEST_H


#include "CEpetra_MultiVector.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Functions AztecOO_StatusTest_Cast() and AztecOO_StatusTest_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the AztecOO_StatusTest table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/*! @name Explicit casting methods */
/*@{*/

CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t AztecOO_StatusTest_Abstract ( 
  CT_AztecOO_StatusTest_ID_t id );

/*@}*/

/*! @name AztecOO_StatusTest destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual AztecOO_StatusTest::~AztecOO_StatusTest()
*/
void AztecOO_StatusTest_Destroy ( 
  CT_AztecOO_StatusTest_ID_t * selfID );

/*@}*/

/*! @name AztecOO_StatusTest member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual bool AztecOO_StatusTest::ResidualVectorRequired() const = 0
*/
boolean AztecOO_StatusTest_ResidualVectorRequired ( 
  CT_AztecOO_StatusTest_ID_t selfID );

/*! @brief Wrapper for 
   virtual AztecOO_StatusType AztecOO_StatusTest::CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated) = 0
*/
CT_AztecOO_StatusType_E_t AztecOO_StatusTest_CheckStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated );

/*! @brief Wrapper for 
   virtual AztecOO_StatusType AztecOO_StatusTest::GetStatus() const = 0
*/
CT_AztecOO_StatusType_E_t AztecOO_StatusTest_GetStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CAZTECOO_STATUSTEST_H */

#endif /* HAVE_CTRILINOS_AZTECOO */


