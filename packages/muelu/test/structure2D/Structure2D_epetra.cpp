/*
 * Structure2D_epetra.cpp
 *
 *  Created on: Oct 24, 2011
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
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"

//
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

/*!
 *  2d structural mechanics example for Epetra
 *
 *  (Nearly) Symmetric problem (except of Dirichlet boundaries) solved with AMG solver using a
 *  3 level multigrid with smoothed aggregation transfer operators.
 *
 *  TODO: There's no matrix amalgamation implemented yet! The example has 1701 nodes, but the graph
 *  of the matrix has 3402 nodes!
 */


int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);


#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  // custom parameters
  LO maxLevels = 2;

  GO maxCoarseSize=1; //FIXME clp doesn't like long long int
  std::string aggOrdering = "natural";
  int minPerAgg=3;
  int maxNbrAlreadySelected=0;

  int globalNumDofs = 3402;
  int nProcs = comm->getSize();
  int nDofsPerNode = 2;

  int nLocalDofs = (int) globalNumDofs / nProcs;
  nLocalDofs = nLocalDofs - (nLocalDofs % nDofsPerNode);
  int nCumulatedDofs = 0;
  sumAll(comm,nLocalDofs, nCumulatedDofs);
  //Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM, 1, nLocalDofs, &nCumulatedDofs );

  if(comm->getRank() == nProcs-1) {
    nLocalDofs += globalNumDofs - nCumulatedDofs;
  }

  std::cout << "PROC: " << comm->getRank() << " numLocalDofs=" << nLocalDofs << std::endl;

  // read in problem
  Epetra_Map emap (globalNumDofs, nLocalDofs, 0, *Xpetra::toEpetra(comm));
  //Epetra_Map emap(3402,0,*Xpetra::toEpetra(comm));
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

  // Epetra_CrsMatrix -> Xpetra::Operator
  RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<CrsOperator> crsOp = Teuchos::rcp(new CrsOperator(exA));
  RCP<Operator> Op = Teuchos::rcp_dynamic_cast<Operator>(crsOp);

  // Epetra_Vector -> Xpetra::Vector
  RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  RCP<MultiVector> xNS = Teuchos::rcp(new Xpetra::EpetraMultiVector(epNS));



  // Epetra_Map -> Xpetra::Map
  const RCP< const Map> map = Xpetra::toXpetra(emap);

  RCP<Hierarchy> H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize(maxCoarseSize);

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);
  Finest->Set("Nullspace",xNS);

  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  dropFact->SetVerbLevel(MueLu::Extreme);
  dropFact->SetFixedBlockSize(nDofsPerNode);
  //RCP<PreDropFunctionConstVal> predrop = rcp(new PreDropFunctionConstVal(0.00001));
  //dropFact->SetPreDropFunction(predrop);
  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory(dropFact));
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
    *out << "aggregate ordering :                    NATURAL" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
    *out << "aggregate ordering :                    RANDOM" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
    *out << "aggregate ordering :                    GRAPH" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  // build transfer operators
  RCP<NullspaceFactory> nspFact = rcp(new NullspaceFactory()); // make sure that we can keep nullspace!!!
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact,nspFact));
  //RCP<PgPFactory> Pfact = rcp( new PgPFactory(TentPFact) );
  //RCP<RFactory> Rfact  = rcp( new GenericRFactory(Pfact));
  RCP<SaPFactory> Pfact  = rcp( new SaPFactory(TentPFact) );
  RCP<RFactory>   Rfact  = rcp( new TransPFactory(Pfact) );
  RCP<RAPFactory> Acfact = rcp( new RAPFactory(Pfact, Rfact) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Finest->Keep("Aggregates",UCAggFact.get());
  Finest->Keep("Nullspace",nspFact.get());

  // build level smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 0.9); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Gauss-Seidel");

  smooProto = Teuchos::rcp( new TrilinosSmoother(Xpetra::UseEpetra,ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  Teuchos::ParameterList status;
  status = H->FullPopulate(*Pfact,*Rfact,*Acfact,*SmooFact,0,maxLevels);

  H->SetCoarsestSolver(*SmooFact,MueLu::PRE);


  *out << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  Finest->print(*out);

  RCP<Level> coarseLevel = H->GetLevel(1);
  coarseLevel->print(*out);

  /*RCP<Level> coarseLevel2 = H->GetLevel(2);
  coarseLevel2->print(*out);
  RCP<MultiVector> nsp2 = coarseLevel2->Get<RCP<MultiVector> >("Nullspace",nspFact.get());
  nsp2->describe(*out,Teuchos::VERB_EXTREME);

  RCP<Level> coarseLevel3 = H->GetLevel(3);
  coarseLevel3->print(*out);*/


  RCP<MultiVector> xLsg = MultiVectorFactory::Build(map,1);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    H->Iterate(*xRhs,10,*xLsg);

    //xLsg->describe(*out,Teuchos::VERB_EXTREME);
  }

  return EXIT_SUCCESS;
}
