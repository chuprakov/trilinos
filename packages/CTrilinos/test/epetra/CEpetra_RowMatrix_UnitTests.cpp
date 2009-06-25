#include "CEpetra_Import.h"
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

/**********************************************************************
void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

