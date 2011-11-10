#include <iostream>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

#include <MueLu_TrilinosSmoother.hpp> //TODO: remove

// Header files defining default types for template parameters.
// These headers must be included after other MueLu/Xpetra headers.
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=int, GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>    // => typedef MueLu::FooClass<Scalar, LocalOrdinal, ...> Foo

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

int main(int argc, char *argv[]) {

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp; //

  //
  // MPI initialization using Teuchos
  //

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  //
  // Parameters
  //

  GlobalOrdinal numGlobalElements = 256; // problem size

  //
  // Construct the problem
  //

  // Construct a Map that puts approximately the same number of equations on each processor
  const Epetra_Map map(numGlobalElements, 0, comm);

  // Get update list and number of local equations from newly created map.
  const size_t         numMyElements    = map.NumMyElements();
  const GlobalOrdinal* myGlobalElements = map.MyGlobalElements();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, map, 3));

  // Add rows one-at-a-time
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {

      //TODO: should be rewritten in an Epetra style
      A->InsertGlobalValues(myGlobalElements[i], 2, 
                            Teuchos::tuple<Scalar> (2.0, -1.0).getRawPtr(),
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] +1).getRawPtr());

    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->InsertGlobalValues(myGlobalElements[i], 2,
                            Teuchos::tuple<Scalar> (-1.0, 2.0).getRawPtr(),
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]).getRawPtr());
    }
    else {
      A->InsertGlobalValues(myGlobalElements[i], 3,
                            Teuchos::tuple<Scalar> (-1.0, 2.0, -1.0).getRawPtr(),
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1).getRawPtr());
    }
  }

  // Complete the fill, ask that storage be reallocated and optimized
  A->FillComplete();

  //
  // Construct a multigrid preconditioner
  //

  // Turns a Epetra_CrsMatrix into a MueLu::Operator
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO, LMO> > mueluA_ = rcp(new Xpetra::EpetraCrsMatrix(A)); //TODO: should not be needed
  RCP<Xpetra::Operator <SC, LO, GO, NO, LMO> > mueluA  = rcp(new Xpetra::CrsOperator<SC, LO, GO, NO, LMO>(mueluA_));

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(mueluA));
  H->setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  FactoryManager M;                         // -
  M.SetFactory("A", rcp(new RAPFactory())); // TODO: to be remove, but will require some work

  H->Setup(M); //Should be instead: H->Setup();

  //
  // Define RHS / LHS
  //

  RCP<Epetra_Vector> X = rcp(new Epetra_Vector(map));
  RCP<Epetra_Vector> B = rcp(new Epetra_Vector(map));
  
  X->PutScalar((Scalar) 0.0);
  B->SetSeed(846930886); B->Random();

#ifndef HAVE_MUELU_BELOS

  //
  // Use AMG directly as an iterative solver (not as a preconditionner)
  //

  int nIts = 9;

  RCP<MueLu::Vector<LO, GO, NO> > mueluX = rcp(new MueLu::Vector<LOG, GO, NO>(X));
  RCP<MueLu::Vector<LO, GO, NO> > mueluB = rcp(new MueLu::Vector<LOG, GO, NO>(B));

  H->Iterate(*mueluB, nIts, *mueluX);

  // Print relative residual norm
  ST::magnitudeType residualNorms = Utils::ResidualNorm(*mueluA, *mueluX, *mueluB)[0];
  if (comm.MyPID() == 0)
    std::cout << "||Residual|| = " << residualNorms << std::endl;

#else // HAVE_MUELU_BELOS

  //
  // Solve Ax = b using AMG as a preconditioner in Belos
  //

  // Operator and Multivector type that will be used with Belos
  typedef Epetra_MultiVector   MV;
  typedef Belos::OperatorT<MV> OP;

  // Define Operator and Preconditioner
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(mueluA)); // Turns a Xpetra::Operator object into a Belos operator
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<SC, LO, GO, NO, LMO>(H));       // Turns a MueLu::Hierarchy object into a Belos operator

  // Construct a Belos LinearProblem object
  RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
  belosProblem->setLeftPrec(belosPrec);
    
  bool set = belosProblem->setProblem();
  if (set == false) {
    std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return EXIT_FAILURE;
  }
    
  // Belos parameter list
  int maxIts = 20;
  double tol = 1e-4;
  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);

  // Create an iterative solver manager
  RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
    
  // Perform solve
  Belos::ReturnType ret = solver->solve();
  
  // Get the number of iterations for this solve.
  std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
  
  // Compute actual residuals.
  int numrhs=1;
  bool badRes = false;
  std::vector<SC> actual_resids(numrhs);
  std::vector<SC> rhs_norm(numrhs);
  RCP<Epetra_MultiVector > resid = rcp(new Epetra_MultiVector(map, numrhs));

  typedef Belos::OperatorTraits<SC, MV, OP> OPT;
  typedef Belos::MultiVecTraits<SC, MV>     MVT;
    
  OPT::Apply(*belosOp, *X, *resid);
  MVT::MvAddMv(-1.0, *resid, 1.0, *B, *resid);
  MVT::MvNorm(*resid, actual_resids);
  MVT::MvNorm(*B, rhs_norm);
  std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  for (int i = 0; i < numrhs; i++) {
    SC actRes = actual_resids[i]/rhs_norm[i];
    std::cout <<"Problem " << i << " : \t" << actRes <<std::endl;
    if (actRes > tol) { badRes = true; }
  }

  // Check convergence
  if (ret != Belos::Converged || badRes) {
    std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;

#endif // HAVE_MUELU_BELOS


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
