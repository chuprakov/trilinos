#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_JADMATRIX_H
#define CEPETRA_JADMATRIX_H


#include "CEpetra_RowMatrix.h"
#include "CEpetra_MultiVector.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_JadMatrix_Cast() and Epetra_JadMatrix_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_JadMatrix table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_JadMatrix_Abstract ( 
  CT_Epetra_JadMatrix_ID_t id );

/* Original C++ prototype:
   Epetra_JadMatrix(const Epetra_RowMatrix & Matrix);
*/
CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID );

/* Original C++ prototype:
   virtual ~Epetra_JadMatrix();
*/
void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID );

/* Original C++ prototype:
   int UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false);
*/
int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure );

/* Original C++ prototype:
   int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;
*/
int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/* Original C++ prototype:
   int ExtractMyEntryView(int CurEntry, double * &Value, int & RowIndex, int & ColIndex);
*/
int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex );

/* Original C++ prototype:
   int ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const;
*/
int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex );

/* Original C++ prototype:
   int NumMyRowEntries(int MyRow, int & NumEntries) const;
*/
int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries );

/* Original C++ prototype:
   int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/* Original C++ prototype:
   int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
*/
int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_JADMATRIX_H */

