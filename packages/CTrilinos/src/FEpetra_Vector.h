
/* True C header file! */


#ifndef FEPETRA_VECTOR_H
#define FEPETRA_VECTOR_H


#include "FEpetra_Map.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef int VectorID;

VectorID Epetra_Vector_Create( MapID mapID );

void Epetra_Vector_Destroy( VectorID vectorID );

void Epetra_Vector_PutScalar( VectorID vectorID, double scalarConstant );

void Epetra_Vector_Random( VectorID vectorID );

void Epetra_Vector_Update(
  VectorID vectorID, double alpha, VectorID vector2ID, double beta );

double Epetra_Vector_Norm2( VectorID vectorID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* FEPETRA_VECTOR_H */
