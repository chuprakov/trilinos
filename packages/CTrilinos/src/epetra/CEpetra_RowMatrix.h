#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_ROWMATRIX_H
#define CEPETRA_ROWMATRIX_H


#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Map.h"
#include "CEpetra_Import.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_RowMatrix_Abstract ( 
  CT_Epetra_RowMatrix_ID_t id );

/* Original C++ prototype:
   virtual ~Epetra_RowMatrix();
*/
void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID );

/* Original C++ prototype:
   virtual int NumMyRowEntries(int MyRow, int & NumEntries) const = 0;
*/
int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries );

/* Original C++ prototype:
   virtual int MaxNumEntries() const = 0;
*/
int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const = 0;
*/
int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/* Original C++ prototype:
   virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const = 0;
*/
int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );

/* Original C++ prototype:
   virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
*/
int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
*/
int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   virtual int InvRowSums(Epetra_Vector& x) const = 0;
*/
int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   virtual int LeftScale(const Epetra_Vector& x) = 0;
*/
int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   virtual int InvColSums(Epetra_Vector& x) const = 0;
*/
int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   virtual int RightScale(const Epetra_Vector& x) = 0;
*/
int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/* Original C++ prototype:
   virtual bool Filled() const = 0;
*/
boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual double NormInf() const = 0;
*/
double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual double NormOne() const = 0;
*/
double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumGlobalNonzeros() const = 0;
*/
int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumGlobalRows() const = 0;
*/
int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumGlobalCols() const= 0;
*/
int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumGlobalDiagonals() const = 0;
*/
int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumMyNonzeros() const = 0;
*/
int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumMyRows() const = 0;
*/
int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumMyCols() const = 0;
*/
int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual int NumMyDiagonals() const = 0;
*/
int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual bool LowerTriangular() const = 0;
*/
boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual bool UpperTriangular() const = 0;
*/
boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Map & RowMatrixRowMap() const = 0;
*/
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Map & RowMatrixColMap() const = 0;
*/
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/* Original C++ prototype:
   virtual const Epetra_Import * RowMatrixImporter() const = 0;
*/
CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_ROWMATRIX_H */

