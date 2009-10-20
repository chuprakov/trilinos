
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


/*! @file CEpetra_Object.h
 * @brief Wrappers for Epetra_Object */

/* True C header file! */


#ifndef CEPETRA_OBJECT_H
#define CEPETRA_OBJECT_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Functions Epetra_Object_Cast() and Epetra_Object_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Object table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/*! @name Explicit casting methods */
/*@{*/

CT_Epetra_Object_ID_t Epetra_Object_Cast ( CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Object_Abstract ( 
  CT_Epetra_Object_ID_t id );

/*@}*/

/*! @name Epetra_Object constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Object::Epetra_Object(int TracebackModeIn = -1, bool set_label = true)
*/
CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );

/*! @brief Wrapper for 
   Epetra_Object::Epetra_Object(const char * const Label, int TracebackModeIn = -1)
*/
CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );

/*! @brief Wrapper for 
   Epetra_Object::Epetra_Object(const Epetra_Object& Object)
*/
CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );

/*@}*/

/*! @name Epetra_Object destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Object::~Epetra_Object()
*/
void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );

/*@}*/

/*! @name Epetra_Object member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual void Epetra_Object::SetLabel(const char * const Label)
*/
void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );

/*! @brief Wrapper for 
   virtual const char * Epetra_Object::Label() const
*/
const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Object::ReportError(const string Message, int ErrorCode) const
*/
int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char Message[], int ErrorCode );

/*@}*/

/*! @name Epetra_Object static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   static void Epetra_Object::SetTracebackMode(int TracebackModeValue)
*/
void Epetra_Object_SetTracebackMode ( int TracebackModeValue );

/*! @brief Wrapper for 
   static int Epetra_Object::GetTracebackMode()
*/
int Epetra_Object_GetTracebackMode (  );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OBJECT_H */

