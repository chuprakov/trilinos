
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


/*! @file CEpetra_Flops.h
 * @brief Wrappers for Epetra_Flops */

/* True C header file! */


#ifndef CEPETRA_FLOPS_H
#define CEPETRA_FLOPS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Functions Epetra_Flops_Cast() and Epetra_Flops_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Flops table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/*! @name Explicit casting methods */
/*@{*/

CT_Epetra_Flops_ID_t Epetra_Flops_Cast ( CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Flops_Abstract ( 
  CT_Epetra_Flops_ID_t id );

/*@}*/

/*! @name Epetra_Flops constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Flops::Epetra_Flops(void)
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );

/*! @brief Wrapper for 
   Epetra_Flops::Epetra_Flops(const Epetra_Flops& Flops_in)
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID );

/*@}*/

/*! @name Epetra_Flops destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Flops::~Epetra_Flops(void)
*/
void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );

/*@}*/

/*! @name Epetra_Flops member wrappers */
/*@{*/

/*! @brief Wrapper for 
   double Epetra_Flops::Flops() const
*/
double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_Flops::ResetFlops()
*/
void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );

/*@}*/

/*! @name Epetra_Flops operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Flops& Epetra_Flops::operator=(const Epetra_Flops& src)
*/
void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_FLOPS_H */

