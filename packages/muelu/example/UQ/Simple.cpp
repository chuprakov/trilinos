#include <iostream>

#include "linear2d_diffusion_scalar_types.hpp"

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"

// Header files defining default types for template parameters.
// These headers must be included after other MueLu/Xpetra headers.
//#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=int, GlobalOrdinal=int
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
#include <MueLu_UseShortNames.hpp>    // => typedef MueLu::FooClass<Scalar, LocalOrdinal, ...> Foo

template <class ST>
class mySimpleClass {
  public: 

  mySimpleClass(ST const & val) : val_(val) {};
  virtual ~mySimpleClass() {};

  private:
    ST val_;

}; //class mySimpleClass

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers

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
/*
#if   defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra; 
#elif defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)  && defined(HAVE_MUELU_AMESOS)
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
#endif
*/
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra; 

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

  // TODO Amesos2 chokes on PCE scalars
  Teuchos::ParameterList coarseParamList;
  coarseParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  coarseParamList.set("relaxation: sweeps", (LO) 1);
  coarseParamList.set("relaxation: damping factor", (double) 1.0);
  RCP<SmootherPrototype> coarsePrototype     = rcp( new TrilinosSmoother("RELAXATION", coarseParamList) );
  RCP<SmootherFactory>   coarseSolverFact      = rcp( new SmootherFactory(coarsePrototype, Teuchos::null) );

  FactoryManager M;
  M.SetFactory("CoarseSolver", coarseSolverFact);

  // Multigrid Hierarchy
  Hierarchy H(A);
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  H.Setup(M);

  //
  // Solve Ax = b
  //

  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);
  
  X->putScalar((Scalar) 0.0);
  B->setSeed(846930886); B->randomize();

  // Use AMG directly as an iterative solver (not as a preconditionner)
  int nIts = 9;

  H.Iterate(*B, nIts, *X);

  // Print relative residual norm
  ST::magnitudeType residualNorms = Utils::ResidualNorm(*A, *X, *B)[0];
  if (comm->getRank() == 0)
    std::cout << "||Residual|| = " << residualNorms << std::endl;

  return EXIT_SUCCESS;
}
