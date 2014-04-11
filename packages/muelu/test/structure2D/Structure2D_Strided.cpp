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
/*
 * Structure2D_visualization.cpp
 *
 *  Created on: Feb 14, 2012
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_CoarseMapFactory.hpp"

//#include "MueLu_AggregationExportFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <MueLu_EpetraOperator.hpp>

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

/*!
 *  2d structural mechanics example for Epetra
 *
 *  (Nearly) Symmetric problem (except of Dirichlet boundaries) solved with AMG solver using a
 *  3 level multigrid with smoothed aggregation transfer operators.
 */

Teuchos::RCP<Vector> runExample(std::vector<size_t> stridingInfo, LocalOrdinal stridedBlockId, GlobalOrdinal offset) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor Mt(myTime);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  // custom parameters
  LO maxLevels = 4;

  GO maxCoarseSize=1; //FIXME clp doesn't like long long int
  std::string aggOrdering = "natural";
  int minPerAgg=3;
  int maxNbrAlreadySelected=0;

  ////////////////////////////////////////////////////////////////////////////////////////
  // prepare redistribution of matrix (parallelization)
  int globalNumDofs = 7020;
  int nProcs = comm->getSize();
  int nDofsPerNode = 2;

  int nLocalDofs = (int) globalNumDofs / nProcs;
  nLocalDofs = nLocalDofs - (nLocalDofs % nDofsPerNode);
  int nCumulatedDofs = 0;
  sumAll(comm,nLocalDofs, nCumulatedDofs);

  if(comm->getRank() == nProcs-1) {
    nLocalDofs += globalNumDofs - nCumulatedDofs;
  }

  std::cout << "PROC: " << comm->getRank() << " numLocalDofs=" << nLocalDofs << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////
  // read in problem
  Epetra_Map emap (globalNumDofs, nLocalDofs, 0, *Xpetra::toEpetra(comm));
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  std::cout << "Reading matrix market file" << std::endl;
  EpetraExt::MatrixMarketFileToCrsMatrix("stru2d_A.txt",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("stru2d_b.txt",emap,ptrf);
  EpetraExt::MatrixMarketFileToMultiVector( "stru2d_ns.txt", emap, ptrNS);
  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

  ////////////////////////////////////////////
  // Epetra_CrsMatrix -> Xpetra::Matrix
  RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
  RCP<Matrix> Op = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);
  Op->SetFixedBlockSize(nDofsPerNode);

  // Epetra_Vector -> Xpetra::Vector
  RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  RCP<MultiVector> xNS = Teuchos::rcp(new Xpetra::EpetraMultiVector(epNS));

  // Epetra_Map -> Xpetra::Map
  const RCP< const Map> map = Xpetra::toXpetra(emap);

  ////////////////////////////////////////////
  // create new MueLu hierarchy
  RCP<Hierarchy> H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize(maxCoarseSize);

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);
  Finest->Set("Nullspace",xNS);

  // prepare CoalesceDropFactory
  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  //dropFact->SetVariableBlockSize();

  // prepare aggregation strategy
  RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
  CoupledAggFact->SetFactory("Graph", dropFact);
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  CoupledAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  CoupledAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
    *out << "aggregate ordering :                    NATURAL" << std::endl;
    CoupledAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
    *out << "aggregate ordering :                    RANDOM" << std::endl;
    CoupledAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
    *out << "aggregate ordering :                    GRAPH" << std::endl;
    CoupledAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  CoupledAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  // build transfer operators
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory());

  RCP<CoarseMapFactory> coarseMapFact = Teuchos::rcp(new CoarseMapFactory());
  coarseMapFact->setStridingData(stridingInfo);
  coarseMapFact->setStridedBlockId(stridedBlockId);
  std::stringstream strOffset;
  strOffset << "{" << offset << "," << offset << "," << offset << "," << offset << "}";
  coarseMapFact->SetParameter("Domain GID offsets",Teuchos::ParameterEntry(strOffset.str()));

  RCP<SaPFactory> Pfact  = rcp( new SaPFactory() );
  //RCP<PgPFactory> Pfact  = rcp( new PgPFactory() );
  //RCP<TentativePFactory> Pfact  = rcp( new TentativePFactory() );
  RCP<Factory>   Rfact  = rcp( new TransPFactory() );
  //RCP<Factory>   Rfact  = rcp( new GenericRFactory() );

  // RAP Factory
  RCP<RAPFactory> Acfact = rcp( new RAPFactory() );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  // register aggregation export factory in RAPFactory
  //RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> > aggExpFact = rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>());
  //aggExpFact->SetParameter("Output filename","aggs_level%LEVELID_proc%PROCID.out");
  //Acfact->AddTransferFactory(aggExpFact);

  // build level smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");

  smooProto = Teuchos::rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  // create coarsest smoother
  RCP<SmootherPrototype> coarsestSmooProto;
  std::string type = "";
  Teuchos::ParameterList coarsestSmooList;
#if defined(HAVE_AMESOS_SUPERLU)
  coarsestSmooProto = Teuchos::rcp( new DirectSolver("Superlu", coarsestSmooList) );
#else
  coarsestSmooProto = Teuchos::rcp( new DirectSolver("Klu", coarsestSmooList) );
#endif
  RCP<SmootherFactory> coarsestSmooFact;
  coarsestSmooFact = rcp(new SmootherFactory(coarsestSmooProto, Teuchos::null));

  FactoryManager M;
  //M.SetFactory("Graph", dropFact);
  //M.SetFactory("UnAmalgamationInfo", dropFact);
  M.SetFactory("Aggregates", CoupledAggFact);
  M.SetFactory("Ptent", TentPFact);
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarsestSmooFact);
  M.SetFactory("CoarseMap", coarseMapFact);

  H->Setup(M, 0, maxLevels);

  RCP<Vector> xLsg = VectorFactory::Build(map);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    H->Iterate(*xRhs,10,*xLsg);

    //xLsg->describe(*out,Teuchos::VERB_EXTREME);
  }

  //
  // Solve Ax = b using AMG as a preconditioner in AztecOO
  //
  {
    RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
    X->PutScalar(0.0);
    Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

    AztecOO aztecSolver(epetraProblem);
    aztecSolver.SetAztecOption(AZ_solver, AZ_cg);

    MueLu::EpetraOperator aztecPrec(H);
    aztecSolver.SetPrecOperator(&aztecPrec);

    int maxIts = 50;
    double tol = 1e-8;

    aztecSolver.Iterate(maxIts, tol);
  }

  return xLsg;
}

}

int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using namespace MueLuTests;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(3);

  // no striding, just one block of size 3, no offset
  Teuchos::RCP<Vector> ref = runExample(stridingInfo, -1, 0);

  int cnt_errors = 0;

  stridingInfo.push_back(1);
  // striding (3,1), use block 0, no offset
  Teuchos::RCP<Vector> lsg2 = runExample(stridingInfo, 0, 0);
  lsg2->update(-1.0, *ref, 1.0);
  if(lsg2->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  stridingInfo.push_back(4);
  // striding (3,1,4), use block 0, no offset
  Teuchos::RCP<Vector> lsg3 = runExample(stridingInfo, 0, 0);
  lsg3->update(-1.0, *ref, 1.0);
  if(lsg3->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  stridingInfo.push_back(3);
  // striding (3,1,4,3), use block 3, no offset
  Teuchos::RCP<Vector> lsg4 = runExample(stridingInfo, 3, 0);
  lsg4->update(-1.0, *ref, 1.0);
  if(lsg4->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  ////////////////////////////////////////////// test with offset
  stridingInfo.clear();
  stridingInfo.push_back(3);

  // no striding, just one block of size 3, offset 35
  Teuchos::RCP<Vector> lsg5 = runExample(stridingInfo, -1, 35);
  lsg5->update(-1.0, *ref, 1.0);
  if(lsg5->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  // no striding, use block 0, offset 35
  Teuchos::RCP<Vector> lsg6 = runExample(stridingInfo, 0, 35);
  lsg6->update(-1.0, *ref, 1.0);
  if(lsg6->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  stridingInfo.push_back(1);
  //striding (3,1), use block 0, offset 35
  Teuchos::RCP<Vector> lsg7 = runExample(stridingInfo, 0, 35);
  lsg7->update(-1.0, *ref, 1.0);
  if(lsg7->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  stridingInfo.push_back(3);
  //striding (3,1,3), use block 2, offset 35
  Teuchos::RCP<Vector> lsg8 = runExample(stridingInfo, 2, 35);
  lsg8->update(-1.0, *ref, 1.0);
  if(lsg8->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  //striding (3,1,3), use block 2, offset 36
  Teuchos::RCP<Vector> lsg9 = runExample(stridingInfo, 2, 36);
  lsg9->update(-1.0, *ref, 1.0);
  if(lsg9->norm2() != Teuchos::ScalarTraits< Scalar >::zero()) { cnt_errors++; }

  if(cnt_errors>0) {
    std::cout << "results do not match. Error" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


