#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "CEpetra_SerialComm.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"

int main( int argc, char* argv[] )
{
  
  /*
   * Data declarations (how old-school is this!)
   */

  int numGlobalElements, numGlobalElements_rtn,numMyElements;

  CT_Epetra_SerialComm_ID_t scommID;
  CT_Epetra_Comm_ID_t commID;

  CT_Epetra_Map_ID_t mapID;
  CT_Epetra_BlockMap_ID_t bmapID;

  CT_Epetra_Vector_ID_t xID, bID;
  CT_Epetra_MultiVector_ID_t mxID, mbID;

  double bnorm, xnorm, expected_bnorm, expected_xnorm, bnorm_err, xnorm_err, err_tol;
  int success = 1;

  /*
   * Executable code
   */
  
  /* Create a serial comm and cast to comm */
  scommID = Epetra_SerialComm_Create();
  commID = Epetra_Comm_Cast(scommID);

  /* Create a map and cast to block map */
  numGlobalElements = 4;
  mapID = Epetra_Map_Create(numGlobalElements, 0, commID);
  bmapID = Epetra_BlockMap_Cast(mapID);

  /* Check the properties of the map */
  numGlobalElements_rtn = Epetra_BlockMap_NumGlobalElements(bmapID);
  printf( "NumGlobalElements = %d\n", numGlobalElements_rtn );
  assert( numGlobalElements == numGlobalElements_rtn );

  numMyElements = Epetra_BlockMap_NumMyElements(bmapID);
  printf( "NumMyElements = %d\n", numMyElements);
  
  /* Create vectors and cast to multivectors */
  xID = Epetra_Vector_Create(bmapID, TRUE);
  bID = Epetra_Vector_Create(bmapID, FALSE);

  mbID = Epetra_MultiVector_Cast(bID);
  mxID = Epetra_MultiVector_Cast(xID);

  /* Do some vector operations */
  Epetra_MultiVector_PutScalar(mbID, 2.0);
  Epetra_MultiVector_Update_WithA(mxID, 2.0, mbID, 0.0); /* x = 2*b */

  Epetra_MultiVector_Norm2(mbID, &bnorm);
  Epetra_MultiVector_Norm2(mxID, &xnorm);

  printf( "2 norm of x = %f\n", xnorm );
  printf( "2 norm of b = %f\n", bnorm );

  /* Test the expected value */
  err_tol = 0.0;
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements );
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements );
  bnorm_err = fabs( expected_bnorm - bnorm ) / expected_bnorm;
  xnorm_err = fabs( expected_xnorm - xnorm ) / expected_xnorm;
  printf( "error in 2 norm of x = %f\n", bnorm_err );
  printf( "error in 2 norm of b = %f\n", xnorm_err );

  if (bnorm_err > err_tol) success = 0;
  if (xnorm_err > err_tol) success = 0;

  /* Clean up memory (in reverse order)! */
  Epetra_MultiVector_Destroy(&mxID);
  Epetra_MultiVector_Destroy(&mbID);
  Epetra_Vector_Destroy(&bID);
  Epetra_Vector_Destroy(&xID);

  Epetra_BlockMap_Destroy(&bmapID);
  Epetra_Map_Destroy(&mapID);

  Epetra_Comm_Destroy(&commID);
  Epetra_SerialComm_Destroy(&scommID);

  /* This should throw an exception and print an error message! */
  /* Epetra_BlockMap_NumGlobalElements(bmapID); */

  if (success)
    printf( "\nEnd Result: TEST PASSED\n" );
  else
    printf( "\nEnd Result: TEST FAILED\n" );
  
  return ( success ? 0 : 1 );

}
