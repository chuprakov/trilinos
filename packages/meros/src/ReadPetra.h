/*
// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef _READPETRA_H_
#define _READPETRA_H_



#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AztecOO.h"

extern int ReadPetraMatrix(Epetra_Map *row_map, Epetra_Map *col_map, 
			   Epetra_CrsMatrix **A, char *fname);

extern int ReadPetraVector(Epetra_Vector *A, char *fname);

extern void ReadAztecVbr(int bsize, int Nnz, int Nbrows,
			 AZ_MATRIX **Amat, char *fname);

/* utility for printing similar to the one in aztecoo only without rhs */

extern void AZ_capt_matrix(double val[], int indx[], int bindx[], int rpntr[],
		       int cpntr[], int bpntr[], 
		      int data_org[]);

#endif


