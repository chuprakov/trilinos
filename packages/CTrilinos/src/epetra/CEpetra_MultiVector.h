#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_MULTIVECTOR_H
#define CEPETRA_MULTIVECTOR_H


#include "CEpetra_BlockMap.h"
#include "Epetra_DataAccess.h"
#include "CEpetra_Vector.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Cast(
  CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut );

/* Original C++ prototype:
   Epetra_MultiVector(const Epetra_MultiVector& Source);
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID );

/* Original C++ prototype:
   Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *A, int MyLDA, int NumVectors);
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, double * A, 
  int MyLDA, int NumVectors );

/* Original C++ prototype:
   Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double **ArrayOfPointers, int NumVectors);
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors );

/* Original C++ prototype:
   Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int *Indices, int NumVectors);
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors );

/* Original C++ prototype:
   Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int StartIndex, int NumVectors);
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors );

/* Original C++ prototype:
   virtual ~Epetra_MultiVector();
*/
void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID );

/* Original C++ prototype:
   int ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_ReplaceGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue );

/* Original C++ prototype:
   int ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/* Original C++ prototype:
   int SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_SumIntoGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue );

/* Original C++ prototype:
   int SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/* Original C++ prototype:
   int ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_ReplaceMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue );

/* Original C++ prototype:
   int ReplaceMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_ReplaceMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/* Original C++ prototype:
   int SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_SumIntoMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue );

/* Original C++ prototype:
   int SumIntoMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
*/
int Epetra_MultiVector_SumIntoMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/* Original C++ prototype:
   int PutScalar (double ScalarConstant);
*/
int Epetra_MultiVector_PutScalar ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant );

/* Original C++ prototype:
   int Random();
*/
int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   int ExtractCopy(double *A, int MyLDA) const;
*/
int Epetra_MultiVector_ExtractCopy_Fill2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double * A, int MyLDA );

/* Original C++ prototype:
   int ExtractCopy(double **ArrayOfPointers) const;
*/
int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers );

/* Original C++ prototype:
   int ExtractView(double **A, int *MyLDA) const;
*/
int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA );

/* Original C++ prototype:
   int ExtractView(double ***ArrayOfPointers) const;
*/
int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers );

/* Original C++ prototype:
   int Dot(const Epetra_MultiVector& A, double *Result) const;
*/
int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result );

/* Original C++ prototype:
   int Abs(const Epetra_MultiVector& A);
*/
int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID );

/* Original C++ prototype:
   int Reciprocal(const Epetra_MultiVector& A);
*/
int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID );

/* Original C++ prototype:
   int Scale(double ScalarValue);
*/
int Epetra_MultiVector_Scale_Self ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarValue );

/* Original C++ prototype:
   int Scale(double ScalarA, const Epetra_MultiVector& A);
*/
int Epetra_MultiVector_Scale ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID );

/* Original C++ prototype:
   int Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis);
*/
int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis );

/* Original C++ prototype:
   int Update(double ScalarA, const Epetra_MultiVector& A, double ScalarB, const Epetra_MultiVector& B, double ScalarThis);
*/
int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis );

/* Original C++ prototype:
   int Norm1 (double * Result) const;
*/
int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/* Original C++ prototype:
   int Norm2 (double * Result) const;
*/
int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/* Original C++ prototype:
   int NormInf (double * Result) const;
*/
int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/* Original C++ prototype:
   int NormWeighted (const Epetra_MultiVector& Weights, double * Result) const;
*/
int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result );

/* Original C++ prototype:
   int MinValue (double * Result) const;
*/
int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/* Original C++ prototype:
   int MaxValue (double * Result) const;
*/
int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/* Original C++ prototype:
   int MeanValue (double * Result) const;
*/
int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/* Original C++ prototype:
   int Multiply(char TransA, char TransB, double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, double ScalarThis );
*/
int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis );

/* Original C++ prototype:
   int Multiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, double ScalarThis );
*/
int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis );

/* Original C++ prototype:
   int ReciprocalMultiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, double ScalarThis );
*/
int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis );

/* Original C++ prototype:
   int SetSeed(unsigned int Seed_in);
*/
int Epetra_MultiVector_SetSeed ( 
  CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in );

/* Original C++ prototype:
   unsigned int Seed();
*/
unsigned int Epetra_MultiVector_Seed ( 
  CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   Epetra_MultiVector& operator = (const Epetra_MultiVector& Source);
*/
void Epetra_MultiVector_Assign ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t SourceID );

/* Original C++ prototype:
   double * const & operator [] (int i) const;
*/
double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i );

/* Original C++ prototype:
   const Epetra_Vector * & operator () (int i) const;
*/
CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i );

/* Original C++ prototype:
   int NumVectors() const;
*/
int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   int MyLength() const;
*/
int Epetra_MultiVector_MyLength ( 
  CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   int GlobalLength() const;
*/
int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   int Stride() const;
*/
int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   bool ConstantStride() const;
*/
boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID );

/* Original C++ prototype:
   int ReplaceMap(const Epetra_BlockMap& map);
*/
int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MULTIVECTOR_H */

