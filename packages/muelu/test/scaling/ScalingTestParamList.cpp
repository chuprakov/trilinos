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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_EasyParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& fancyout = *fancy;
  fancyout.setOutputToRootOnly(0);


  // =========================================================================
  // Parameters initialization
  // =========================================================================
  Teuchos::CommandLineProcessor clp(false);

  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName       = "scalingTest.xml"; clp.setOption("xml",                   &xmlFileName,      "read parameters from a file [default = 'scalingTest.xml']");
  bool        printTimings      = true;              clp.setOption("timings", "notimings",  &printTimings,     "print timings to screen");
  int         writeMatricesOPT  = -2;                clp.setOption("write",                 &writeMatricesOPT, "write matrices to file (-1 means all; i>=0 means level i)");
  std::string solveType         = "cg";              clp.setOption("solver",                &solveType,        "solve type: (none | cg | gmres | standalone)");
  double      tol               = 1e-12;             clp.setOption("tol",                   &tol,              "solver convergence tolerance");

  std::string mapFile;                               clp.setOption("map",                   &mapFile,          "map data file");
  std::string matrixFile;                            clp.setOption("matrix",                &matrixFile,       "matrix data file");
  std::string coordFile;                             clp.setOption("coords",                &coordFile,        "coordinates data file");
  int         numRebuilds       = 0;                 clp.setOption("rebuild",               &numRebuilds, "#times to rebuild hierarchy");

  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  // =========================================================================
  // Problem construction
  // =========================================================================
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time"))), tm;

  comm->barrier();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

  RCP<Matrix>      A;
  RCP<MultiVector> coordinates;
  RCP<MultiVector> nullspace;
  if (matrixFile.empty()) {
    fancyout << "========================================================\n" << xpetraParameters << matrixParameters;

    // Retrieve matrix parameters (they may have been changed on the command line), and pass them to Galeri.
    // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    // A perfect distribution is only possible when the #processors is a perfect square.
    // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
    // size. For example, np=14 will give a 7-by-2 distribution.
    // If you don't want Galeri to do this, specify mx or my on the galeriList.
    Teuchos::ParameterList pl = matrixParameters.GetParameterList();
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", pl.get("nx",nx));
    galeriList.set("ny", pl.get("ny",ny));
    galeriList.set("nz", pl.get("nz",nz));
    // galeriList.set("mx", comm->getSize());
    // galeriList.set("my", 1);

    // Create map and coordinates
    // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
    // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
    RCP<const Map> map;
    if (matrixParameters.GetMatrixType() == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",map,matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace2D" || matrixParameters.GetMatrixType() == "Star2D" || matrixParameters.GetMatrixType() == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D",map,matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace3D" || matrixParameters.GetMatrixType() == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D",map,matrixParameters.GetParameterList());
    }
    // Expand map to do multiple DOF per node for block problems
    if (matrixParameters.GetMatrixType() == "Elasticity2D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 2);
    if (matrixParameters.GetMatrixType() == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 3);

    if (comm->getRank() == 0) {
      GO mx = galeriList.get("mx", -1);
      GO my = galeriList.get("my", -1);
      GO mz = galeriList.get("mz", -1);
      fancyout << "Processor subdomains in x direction: " << mx << std::endl
          << "Processor subdomains in y direction: " << my << std::endl
          << "Processor subdomains in z direction: " << mz << std::endl
          << "========================================================" << std::endl;
    }

    Teuchos::ParameterList matrixParams = matrixParameters.GetParameterList();
    matrixParams.set("mx", galeriList.get("mx", -1));
    matrixParams.set("my", galeriList.get("my", -1));
    matrixParams.set("mz", galeriList.get("mz", -1));
    if (matrixParameters.GetMatrixType() == "Elasticity2D" || matrixParameters.GetMatrixType() == "Elasticity3D") {
      // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
      matrixParams.set("right boundary" , "Neumann");
      matrixParams.set("bottom boundary", "Neumann");
      matrixParams.set("top boundary"   , "Neumann");
      matrixParams.set("front boundary" , "Neumann");
      matrixParams.set("back boundary"  , "Neumann");
    }

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParams);
    A = Pr->BuildMatrix();

    nullspace = MultiVectorFactory::Build(map, 1);
    if (matrixParameters.GetMatrixType() == "Elasticity2D" ||
        matrixParameters.GetMatrixType() == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((matrixParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);

    } else {
      nullspace->putScalar(one);
    }
  } else {
    RCP<const Map> map = (mapFile.empty() ? Teuchos::null : Utils2::ReadMap(mapFile, xpetraParameters.GetLib(), comm));
    comm->barrier();

    A = Utils::Read(matrixFile, map);
    comm->barrier();

    coordinates = Utils2::ReadMultiVector(coordFile, map);

    nullspace = MultiVectorFactory::Build(map, 1);
    nullspace->putScalar(one);
  }
  RCP<const Map> map = A->getRowMap();

  comm->barrier();
  tm = Teuchos::null;

  fancyout << "Galeri complete.\n========================================================" << std::endl;

  // =========================================================================
  // Preconditioner construction
  // =========================================================================
  comm->barrier();
  bool useEasy = true;
  {
    // XML files for the original interpreter always contain "Hierarchy" sublist
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);
    if (paramList.isSublist("Hierarchy"))
      useEasy = false;
  }

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));
  RCP<HierarchyManager> mueLuFactory;
  if (useEasy == false)
    mueLuFactory = rcp(new ParameterListInterpreter(xmlFileName, *comm));
  else
    mueLuFactory = rcp(new EasyParameterListInterpreter(xmlFileName, *comm));

  comm->barrier();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

  RCP<Hierarchy> H;
  for (int i = 0; i <= numRebuilds; i++) {
    A->SetMaxEigenvalueEstimate(-one);

    H = mueLuFactory->CreateHierarchy();
    H->GetLevel(0)->Set("A",           A);
    H->GetLevel(0)->Set("Nullspace",   nullspace);
    H->GetLevel(0)->Set("Coordinates", coordinates);
    mueLuFactory->SetupHierarchy(*H);
  }

  comm->barrier();
  tm = Teuchos::null;

  // =========================================================================
  // System solution (Ax = b)
  // =========================================================================
  comm->barrier();
  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));

  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);

  {
    // we set seed for reproducibility
    X->setSeed(846930886);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(1.0/norms[0]);
    X->putScalar(zero);
  }
  tm = Teuchos::null;

  if (writeMatricesOPT > -2)
    H->Write(writeMatricesOPT, writeMatricesOPT);

  comm->barrier();
  if (solveType == "none") {
    // Do not perform a solve

  } else if (solveType == "standalone") {
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Fixed Point Solve")));

    H->IsPreconditioner(false);
    H->Iterate(*B, 25, *X);

  } else if (solveType == "cg" || solveType == "gmres") {
#ifdef HAVE_MUELU_BELOS
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Belos Solve")));

    // Operator and Multivector type that will be used with Belos
    typedef MultiVector          MV;
    typedef Belos::OperatorT<MV> OP;

    H->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(A)); // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp <SC, LO, GO, NO, LMO>(H)); // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    int maxIts = 2000;
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
    belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency",      1);
    belosList.set("Output Style",          Belos::Brief);

    // Create an iterative solver manager
    RCP< Belos::SolverManager<SC, MV, OP> > solver;
    if (solveType == "cg")
      solver = rcp(new Belos::BlockCGSolMgr   <SC, MV, OP>(belosProblem, rcp(&belosList, false)));
    else if (solveType == "gmres")
      solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

    // Perform solve
    Belos::ReturnType ret = Belos::Unconverged;
    try {
      ret = solver->solve();

      // Get the number of iterations for this solve.
      fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

    } catch(...) {
      fancyout << std::endl << "ERROR:  Belos threw an error! " << std::endl;
    }

    // Check convergence
    if (ret != Belos::Converged)
      fancyout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    else
      fancyout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
#endif //ifdef HAVE_MUELU_BELOS
  } else {
    throw MueLu::Exceptions::RuntimeError("Unknown solver type: \"" + solveType + "\"");
  }
  comm->barrier();
  tm = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  if (printTimings) {
    TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union);
    MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();
  }

  return 0;
} //main
