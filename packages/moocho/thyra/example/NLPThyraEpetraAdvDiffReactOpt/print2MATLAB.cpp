#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_CrsMatrix.h"
#include <iostream>

/* ======== ================ *
 * function CrsMatrix2MATLAB *
 * ======== ================ *
 *
 * Print out a CrsMatrix in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 *
 * Return code:        true if matrix has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to the distributed CrsMatrix to 
 *                     print out
 * - ostream &         reference to output stream 
 */

bool CrsMatrix2MATLAB(const Epetra_CrsMatrix & A, ostream & outfile) 

{

  int MyPID = A.Comm().MyPID(); 
  int NumProc = A.Comm().NumProc();

  // work only on transformed matrices;
  if( A.IndicesAreLocal() == false ) {
    if( MyPID == 0 ) { 
      cerr << "ERROR in "<< __FILE__ << ", line " << __LINE__ << endl;
      cerr << "Function CrsMatrix2MATLAB accepts\n";
      cerr << "transformed matrices ONLY. Please call A.TransformToLoca()\n";
      cerr << "on your matrix A to that purpose.\n";
      cerr << "Now returning...\n";
    }
    return false;
  }

  int NumMyRows = A.NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalNonzeros = A.NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = A.IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  // write on file the dimension of the matrix

  if( MyPID==0 ) {
    outfile << "A = spalloc(";
    outfile << NumGlobalRows << ',' << NumGlobalRows;
    outfile << ',' << NumGlobalNonzeros << ");\n";
  }

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    A.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "\n\n% On proc " << Proc << ": ";
      outfile << NumMyRows << " rows and ";
      outfile << A.NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {

        GlobalRow = A.GRID(MyRow);

        NumNzRow = A.NumMyEntries(MyRow);
        double *Values = new double[NumNzRow];
        int *Indices = new int[NumNzRow];

        A.ExtractMyRowCopy(MyRow, NumNzRow, 
                           NumEntries, Values, Indices);
        // print out the elements with MATLAB syntax
        for( int j=0 ; j<NumEntries ; ++j ) {
          outfile << "A(" << GlobalRow  + IndexBase 
                  << "," << A.GCID(Indices[j]) + IndexBase
                  << ") = " << Values[j] << ";\n";
        }

        delete Values;
        delete Indices;
      }
      
    }
    A.Comm().Barrier();
  }

  return true;

}


/* ======== ============= *
 * function Vector2MATLAB *
 * ======== ============= *
 *
 * Print out a Epetra_Vector in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 * Return code:        true if vector has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_Vector     reference to vector
 * - ostream &         reference to output stream 
 */

bool Vector2MATLAB( const Epetra_Vector & v, ostream & outfile)
{
  
  int MyPID = v.Comm().MyPID(); 
  int NumProc = v.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();
  
  // print out on cout if no filename is provided

  // write on file the dimension of the matrix

  if( MyPID == 0 ) outfile << "v = zeros(" << GlobalLength << ",1)\n";

  int NumMyElements = v.Map().NumMyElements();
  // get update list
  int * MyGlobalElements = v.Map().MyGlobalElements( );
  
  int Row;

  int IndexBase = v.Map().IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    v.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "% On proc " << Proc << ": ";
      outfile << MyLength << " rows of ";
      outfile << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        outfile << "v(" << MyGlobalElements[Row] + IndexBase
             << ") = " << v[Row] << ";\n";
      }
      
    }
      
    v.Comm().Barrier();
  }

  return true;

} /* Vector2MATLAB */


/* ======== =============== *
 * function FEVector2MATLAB *
 * ======== =============== *
 *
 * Print out a Epetra_Vector in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 * Return code:        true if vector has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_FEVector   reference to FE vector
 * - ostream &         reference to output stream 
 */

bool FEVector2MATLAB( const Epetra_FEVector & v, ostream & outfile)
{
  
  int MyPID = v.Comm().MyPID(); 
  int NumProc = v.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();
  
  // print out on cout if no filename is provided

  // write on file the dimension of the matrix

  if( MyPID == 0 ) outfile << "v = zeros(" << GlobalLength << ",1)\n";

  int NumMyElements = v.Map().NumMyElements();
  // get update list
  int * MyGlobalElements = v.Map().MyGlobalElements( );
  
  int Row;

  int IndexBase = v.Map().IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    v.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "% On proc " << Proc << ": ";
      outfile << MyLength << " rows of ";
      outfile << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        outfile << "v(" << MyGlobalElements[Row] + IndexBase
             << ") = " << v[0][Row] << ";\n";
      }
      
    }
      
    v.Comm().Barrier();
  }

  return true;

} /* FEVector2MATLAB */

