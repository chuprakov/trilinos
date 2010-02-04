
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


/*! @file CTrilinos_flex_enums.h
 * @brief Defines structs and enums needed for CTrilinos. */


#ifndef CTRILINOS_FLEX_ENUMS_H
#define CTRILINOS_FLEX_ENUMS_H


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef union {
  CTrilinos_Universal_ID_t universal;
  CT_Epetra_SerialComm_ID_t Epetra_SerialComm;
  CT_Epetra_Comm_ID_t Epetra_Comm;
} CT_Epetra_SerialComm_ID_Flex_t;

typedef union {
  CTrilinos_Universal_ID_t universal;
  CT_Epetra_Map_ID_t Epetra_Map;
  CT_Epetra_BlockMap_ID_t Epetra_BlockMap;
} CT_Epetra_Map_ID_Flex_t;

typedef union {
  CTrilinos_Universal_ID_t universal;
  CT_Epetra_Vector_ID_t Epetra_Vector;
  CT_Epetra_MultiVector_ID_t Epetra_MultiVector;
} CT_Epetra_Vector_ID_Flex_t;


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
