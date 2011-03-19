
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



/*! @file CPliris.h
 * @brief Wrappers for Pliris */

/* True C header file! */


#ifndef CPLIRIS_H
#define CPLIRIS_H

#ifdef HAVE_CTRILINOS_PLIRIS


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Pliris constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Pliris::Pliris(Epetra_Vector * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_Pliris_ID_t Pliris_Create ( 
  CT_Epetra_Vector_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   Pliris::Pliris()
*/
CT_Pliris_ID_t Pliris_Create_Default (  );

/*@}*/

/*! @name Pliris destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Pliris::~Pliris(void)
*/
void Pliris_Destroy ( CT_Pliris_ID_t * selfID );

/*@}*/

/*! @name Pliris member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Pliris::SetLHS(Epetra_MultiVector * X)
*/
int Pliris_SetLHS ( 
  CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t XID );

/*! @brief Wrapper for 
   int Pliris::SetRHS(Epetra_MultiVector * B)
*/
int Pliris_SetRHS ( 
  CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   int Pliris::SetMatrix(Epetra_Vector * A)
*/
int Pliris_SetMatrix ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID );

/*! @brief Wrapper for 
   int Pliris::SetMatrix(Epetra_SerialDenseVector * A)
*/
int Pliris_SetMatrix_Serial ( 
  CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AID );

/*! @brief Wrapper for 
   int Pliris::GetDistribution( int * nprocs_row, int * number_of_unknowns, int * nrhs, int * my_rows, int * my_cols, int * my_first_row, int * my_first_col, int * my_rhs, int * my_row, int * my_col )
*/
int Pliris_GetDistribution ( 
  CT_Pliris_ID_t selfID, int * nprocs_row, int * number_of_unknowns, 
  int * nrhs, int * my_rows, int * my_cols, int * my_first_row, 
  int * my_first_col, int * my_rhs, int * my_row, int * my_col );

/*! @brief Wrapper for 
   int Pliris::FactorSolve( Epetra_Vector * A, int my_rows, int my_cols, int* matrix_size, int* num_procsr, int* num_rhs, double* secs)
*/
int Pliris_FactorSolve ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, int my_rows, 
  int my_cols, int * matrix_size, int * num_procsr, int * num_rhs, 
  double * secs );

/*! @brief Wrapper for 
   int Pliris::FactorSolve( Epetra_SerialDenseVector * AA, int my_rows, int my_cols, int* matrix_size, int* num_procsr, int* num_rhs, double* secs)
*/
int Pliris_FactorSolve_Serial ( 
  CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AAID, 
  int my_rows, int my_cols, int * matrix_size, int * num_procsr, 
  int * num_rhs, double * secs );

/*! @brief Wrapper for 
   int Pliris::Factor( Epetra_Vector* A, int* matrix_size, int* num_procsr, int* permute, double* secs)
*/
int Pliris_Factor ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, 
  int * matrix_size, int * num_procsr, int * permute, 
  double * secs );

/*! @brief Wrapper for 
   int Pliris::Solve(int* permute, int* num_rhs)
*/
int Pliris_Solve ( 
  CT_Pliris_ID_t selfID, int * permute, int * num_rhs );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_PLIRIS */

#endif /* CPLIRIS_H */

