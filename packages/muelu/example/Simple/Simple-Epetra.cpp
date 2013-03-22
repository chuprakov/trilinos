// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
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

  // Turns a Epetra_CrsMatrix into a MueLu::Matrix
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO, LMO> > mueluA_ = rcp(new Xpetra::EpetraCrsMatrix(A)); //TODO: should not be needed
  RCP<Xpetra::Matrix <SC, LO, GO, NO, LMO> > mueluA  = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO, LMO>(mueluA_));

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(mueluA));
  H->setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  H->Setup();

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

  // Wrap Epetra Vectors into Xpetra Vectors
  RCP<Vector> mueluX = rcp(new Xpetra::EpetraVector(X));
  RCP<Vector> mueluB = rcp(new Xpetra::EpetraVector(B));

  H->Iterate(*mueluB, nIts, *mueluX);

  // Print relative residual norm
  ST::magnitudeType residualNorms = Utils::ResidualNorm(*mueluA, *mueluX, *mueluB)[0];
  if (comm.MyPID() == 0)
    std::cout << "||Residual|| = " << residualNorms << std::endl;

#else // HAVE_MUELU_BELOS

  //
  // Solve Ax = b using AMG as a preconditioner in Belos
  //

  // Matrix and Multivector type that will be used with Belos
  typedef Epetra_MultiVector   MV;
  typedef Belos::OperatorT<MV> OP;

  // Define Operator and Preconditioner
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(mueluA)); // Turns a Xpetra::Matrix object into a Belos operator
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
