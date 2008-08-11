#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include <vector>

using namespace Galeri;

bool TestAmesos(char ProblemType[],
		Epetra_RowMatrix& A,
		int NumVectors)
{

  const Epetra_BlockMap& Map = A.OperatorDomainMap();

  Epetra_MultiVector x2(Map,NumVectors);
  Epetra_MultiVector x1(Map,NumVectors);
  Epetra_MultiVector x(Map,NumVectors);
  Epetra_MultiVector b(Map,NumVectors);
  Epetra_MultiVector residual(Map,NumVectors);
  Epetra_MultiVector temp(Map,NumVectors);
  
  Teuchos::ParameterList ParamList ;
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;

  // Note that Abase is created with an empty Problem, none of A, x or b
  // have been specified at this point.  
  Abase = Afactory.Create(ProblemType,Problem ); 
  assert (Abase != 0);

  //  Factor A
  Problem.SetOperator(&A);
  EPETRA_CHK_ERR(Abase->SymbolicFactorization()); 
  EPETRA_CHK_ERR(Abase->NumericFactorization()); 

  b.Random();

  //  Solve Ax = b 
  //
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);
  EPETRA_CHK_ERR(Abase->Solve()); 

  //  Solve Ax1 = x 
  //
  Problem.SetLHS(&x1);
  Problem.SetRHS(&x);
  EPETRA_CHK_ERR(Abase->Solve()); 

  //  Solve Ax2 = x1
  //
  Problem.SetLHS(&x2);
  Problem.SetRHS(&x1);
  EPETRA_CHK_ERR(Abase->Solve()); 

  //  Compute the residual: A^3 x2 - b
  //
  A.Multiply(false,x2, temp) ; //  temp = A x2
  A.Multiply(false,temp, x2) ; //  x2 = A^2 x2
  A.Multiply(false,x2, temp) ; //  temp = A^3 x2
  residual.Update( 1.0, temp, -1.0, b, 0.0 ) ;

  double norm_residual ;
  residual.Norm2( &norm_residual ) ; 

  if (A.Comm().MyPID() == 0) {
    cout << "norm2(A^3 x-b) = " << norm_residual << endl ; 
  }

  delete Abase;

  if (norm_residual < (1e-5))
    return(true);
  else
    return(false);
  
}

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumVectors = 2;

  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", 8);
  GaleriList.set("ny", 8 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  Amesos Factory;  
  
  vector<string> SolverType;
  //  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Pardiso");
  SolverType.push_back("Amesos_Taucs");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");
  SolverType.push_back("Amesos_Dscpack");

  bool TestPassed = true;

  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) {
    string Solver = SolverType[i];
        cout  << Solver << " next  " << endl;
    if (Factory.Query((char*)Solver.c_str())) {
      if (Comm.MyPID() == 0)
        cout << "Testing " << Solver << endl;
      if(TestAmesos((char*)Solver.c_str(), *A, NumVectors) == false) {
        cout  << Solver << " Failed " << endl;
	TestPassed = false;
      } else { 
        cout  << Solver << " Passed " << endl;
      } 
    } else
      if (Comm.MyPID() == 0) {
	cerr << endl;
	cerr << "WARNING: SOLVER `" << Solver << "' NOT TESTED" << endl;
	cerr << endl;
      }
  }
   
  delete A;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (TestPassed) {
    if (Comm.MyPID() == 0) 
      cout << "TESTS PASSED!" << endl;
    return( EXIT_SUCCESS );
  } 
  else {
    if (Comm.MyPID() == 0) 
      cout << "TESTS FAILED!" << endl;
    return( EXIT_FAILURE );
  }

}
