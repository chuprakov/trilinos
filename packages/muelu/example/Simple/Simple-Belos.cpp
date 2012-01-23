//
// To compile and run this example, the optional dependency Belos must be available.
//

#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

#include <MueLu_TrilinosSmoother.hpp> //TODO: remove
#include <MueLu_RAPFactory.hpp> //TMP

// Header files defining default types for template parameters.
// These headers must be included after other MueLu/Xpetra headers.
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=int, GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>    // => typedef MueLu::FooClass<Scalar, LocalOrdinal, ...> Foo

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp; //

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  // problem size
  GlobalOrdinal numGlobalElements = 256;

  // linear algebra library
  // this example require Tpetra+Ifpack2+Amesos2 or Epetra+Ifpack+Amesos
#if   defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra; 
#elif defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)  && defined(HAVE_MUELU_AMESOS)
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
#endif

  //
  // Construct the problem
  //

  // Construct a Map that puts approximately the same number of equations on each processor
  RCP<const Map> map = MapFactory::createUniformContigMap(lib, numGlobalElements, comm);

  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getNodeElementList();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<Operator> A = rcp(new CrsOperator(map, 3));

  // Add rows one-at-a-time
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<Scalar> (2.0, -1.0));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]), 
                            Teuchos::tuple<Scalar> (-1.0, 2.0));
    }
    else {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<Scalar> (-1.0, 2.0, -1.0));
    }
  }

  // Complete the fill, ask that storage be reallocated and optimized
  A->fillComplete();

  //
  // Construct a multigrid preconditioner
  //

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(A));
  H->setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  FactoryManager M;                         // -
  M.SetFactory("A", rcp(new RAPFactory())); // TODO: to be remove, but will require some work

  H->Setup(M); //Should be instead: H->Setup();

  //
  // Solve Ax = b using AMG as a preconditioner in Belos
  //

  // Define RHS / LHS
  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);
  X->putScalar((Scalar) 0.0);
  B->setSeed(846930886); B->randomize();
  //B->putScalar((Scalar) 1.0);

  // Operator and Multivector type that will be used with Belos
  typedef MultiVector          MV;
  typedef Belos::OperatorT<MV> OP;

  // Define Operator and Preconditioner
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(A)); // Turns a Xpetra::Operator object into a Belos operator
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<SC, LO, GO, NO, LMO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator

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
  RCP<MultiVector> resid = MultiVectorFactory::Build(map, numrhs); 

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

  return EXIT_SUCCESS;
}

//TODO: check results
//TODO: rename rhs_norm etc. and simplify computation of residual
//TODO: output in parallel
