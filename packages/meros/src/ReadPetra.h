#ifndef _READPETRA_H_
#define _READPETRA_H_



#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AztecOO.h"

extern int ReadPetraMatrix(Epetra_Map *map, 
			   Epetra_CrsMatrix **A, char *fname);

extern int ReadPetraVector(Epetra_Vector *A, char *fname);

extern void ReadAztecVbr(int bsize, int Nnz, int Nbrows,
			 AZ_MATRIX **Amat, char *fname);

/* utility for printing similar to the one in aztecoo only without rhs */

extern void AZ_capt_matrix(double val[], int indx[], int bindx[], int rpntr[],
		       int cpntr[], int bpntr[], 
		      int data_org[]);

#endif


