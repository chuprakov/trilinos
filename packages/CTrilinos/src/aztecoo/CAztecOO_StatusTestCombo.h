
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



/*! @file CAztecOO_StatusTestCombo.h
 * @brief Wrappers for AztecOO_StatusTestCombo */

/* True C header file! */


#ifndef CAZTECOO_STATUSTESTCOMBO_H
#define CAZTECOO_STATUSTESTCOMBO_H


#include "CAztecOO_StatusTest.h"
#include "CEpetra_MultiVector.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Functions AztecOO_StatusTestCombo_Cast() and AztecOO_StatusTestCombo_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the AztecOO_StatusTestCombo table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/*! @name Explicit casting methods */
/*@{*/

CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t AztecOO_StatusTestCombo_Abstract ( 
  CT_AztecOO_StatusTestCombo_ID_t id );

/*@}*/

/*! @name AztecOO_StatusTestCombo constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t)
*/
CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create ( 
  CT_ComboType_E_t t );

/*! @brief Wrapper for 
   AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a)
*/
CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_OneTest ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID );

/*! @brief Wrapper for 
   AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a, AztecOO_StatusTest& b)
*/
CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_TwoTests ( 
  CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID, 
  CT_AztecOO_StatusTest_ID_t bID );

/*@}*/

/*! @name AztecOO_StatusTestCombo destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual AztecOO_StatusTestCombo::~AztecOO_StatusTestCombo()
*/
void AztecOO_StatusTestCombo_Destroy ( 
  CT_AztecOO_StatusTestCombo_ID_t * selfID );

/*@}*/

/*! @name AztecOO_StatusTestCombo member wrappers */
/*@{*/

/*! @brief Wrapper for 
   AztecOO_StatusTestCombo& AztecOO_StatusTestCombo::AddStatusTest(AztecOO_StatusTest& a)
*/
CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_AddStatusTest ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID, 
  CT_AztecOO_StatusTest_ID_t aID );

/*! @brief Wrapper for 
   bool AztecOO_StatusTestCombo::ResidualVectorRequired() const
*/
boolean AztecOO_StatusTestCombo_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID );

/*! @brief Wrapper for 
   AztecOO_StatusType AztecOO_StatusTestCombo::CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated)
*/
CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_CheckStatus ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated );

/*! @brief Wrapper for 
   AztecOO_StatusType AztecOO_StatusTestCombo::GetStatus() const
*/
CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_GetStatus ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID );

/*! @brief Wrapper for 
   ComboType AztecOO_StatusTestCombo::GetComboType() const
*/
CT_ComboType_E_t AztecOO_StatusTestCombo_GetComboType ( 
  CT_AztecOO_StatusTestCombo_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CAZTECOO_STATUSTESTCOMBO_H */

#endif /* HAVE_CTRILINOS_AZTECOO */


