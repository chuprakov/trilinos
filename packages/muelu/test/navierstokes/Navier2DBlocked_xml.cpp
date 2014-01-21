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
 * Navier2D_epetra.cpp
 *
 *  Created on: Mar 26, 2011
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
#include <Teuchos_ArrayRCP.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
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
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include <MueLu_ParameterListInterpreter.hpp>

#include "MueLu_UseDefaultTypes.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>


namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

// helper routines
  bool SplitMatrix2x2(Teuchos::RCP<const Epetra_CrsMatrix> A,
                      const Epetra_Map& A11rowmap,
                      const Epetra_Map& A22rowmap,
                      Teuchos::RCP<Epetra_CrsMatrix>& A11,
                      Teuchos::RCP<Epetra_CrsMatrix>& A12,
                      Teuchos::RCP<Epetra_CrsMatrix>& A21,
                      Teuchos::RCP<Epetra_CrsMatrix>& A22)
  {
    if (A==Teuchos::null)
      {
        std::cout << "ERROR: SplitMatrix2x2: A==null on entry" << std::endl;
        return false;
      }

    const Epetra_Comm& Comm   = A->Comm();
    const Epetra_Map&  A22map = A22rowmap;
    const Epetra_Map&  A11map = A11rowmap;

    //----------------------------- create a parallel redundant map of A22map
    std::map<int,int> a22gmap;
    {
      std::vector<int> a22global(A22map.NumGlobalElements());
      std::cout << "numGlobal = " << A22map.NumGlobalElements() << std::endl;
      int count=0;
      for (int proc=0; proc<Comm.NumProc(); ++proc)
        {
          int length = 0;
          if (proc==Comm.MyPID())
            {
              std::cout << "numLocal = " << A22map.NumMyElements() << std::endl;
              for (int i=0; i<A22map.NumMyElements(); ++i)
                {
                  a22global[count+length] = A22map.GID(i);
                  ++length;
                }
            }
          Comm.Broadcast(&length,1,proc);
          Comm.Broadcast(&a22global[count],length,proc);
          count += length;
        }
      if (count != A22map.NumGlobalElements())
        {
          std::cout << "ERROR SplitMatrix2x2: mismatch in dimensions" << std::endl;
          return false;
        }

      // create the map
      for (int i=0; i<count; ++i)
        a22gmap[a22global[i]] = 1;
      a22global.clear();
    }

    //--------------------------------------------------- create matrix A22
    A22 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A22map,100));
    {
      std::vector<int>    a22gcindices(100);
      std::vector<double> a22values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A22map.MyGID(grid)==false)
            continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << std::endl;
              return false;
            }

          if (numentries>(int)a22gcindices.size())
            {
              a22gcindices.resize(numentries);
              a22values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid in a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr==a22gmap.end()) continue;
              //std::cout << gcid << " ";
              a22gcindices[count] = gcid;
              a22values[count]    = values[j];
              ++count;
            }
          //std::cout << std::endl; fflush(stdout);
          // add this filtered row to A22
          err = A22->InsertGlobalValues(grid,count,&a22values[0],&a22gcindices[0]);
          if (err<0)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << std::endl;
              return false;
            }

        } //for (int i=0; i<A->NumMyRows(); ++i)
      a22gcindices.clear();
      a22values.clear();
    }
    A22->FillComplete();
    A22->OptimizeStorage();

    //----------------------------------------------------- create matrix A11
    A11 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A11map,100));
    {
      std::vector<int>    a11gcindices(100);
      std::vector<double> a11values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A11map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << std::endl;
              return false;
            }

          if (numentries>(int)a11gcindices.size())
            {
              a11gcindices.resize(numentries);
              a11values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid as part of a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr!=a22gmap.end()) continue;
              a11gcindices[count] = gcid;
              a11values[count] = values[j];
              ++count;
            }
          err = A11->InsertGlobalValues(grid,count,&a11values[0],&a11gcindices[0]);
          if (err<0)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << std::endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a11gcindices.clear();
      a11values.clear();
    }
    A11->FillComplete();
    A11->OptimizeStorage();

    //---------------------------------------------------- create matrix A12
    A12 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A11map,100));
    {
      std::vector<int>    a12gcindices(100);
      std::vector<double> a12values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A11map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << std::endl;
              return false;
            }

          if (numentries>(int)a12gcindices.size())
            {
              a12gcindices.resize(numentries);
              a12values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid as part of a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr==a22gmap.end()) continue;
              a12gcindices[count] = gcid;
              a12values[count] = values[j];
              ++count;
            }
          err = A12->InsertGlobalValues(grid,count,&a12values[0],&a12gcindices[0]);
          if (err<0)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << std::endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a12values.clear();
      a12gcindices.clear();
    }
    A12->FillComplete(A22map,A11map);
    A12->OptimizeStorage();

    //----------------------------------------------------------- create A21
    A21 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A22map,100));
    {
      std::vector<int>    a21gcindices(100);
      std::vector<double> a21values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A22map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << std::endl;
              return false;
            }

          if (numentries>(int)a21gcindices.size())
            {
              a21gcindices.resize(numentries);
              a21values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid as part of a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr!=a22gmap.end()) continue;
              a21gcindices[count] = gcid;
              a21values[count] = values[j];
              ++count;
            }
          err = A21->InsertGlobalValues(grid,count,&a21values[0],&a21gcindices[0]);
          if (err<0)
            {
              std::cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << std::endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a21values.clear();
      a21gcindices.clear();
    }
    A21->FillComplete(A11map,A22map);
    A21->OptimizeStorage();

    //-------------------------------------------------------------- tidy up
    a22gmap.clear();
    return true;
  }

}

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */


int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLuTests;
  using namespace Teuchos;

  typedef Xpetra::StridedMap<int,int>        StridedMap;
  typedef Xpetra::StridedMapFactory<int,int> StridedMapFactory;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  // default parameters
  std::string xmlFile = "myXML.xml";

  // Note: use --help to list available options.
  CommandLineProcessor clp(false);
  clp.setOption("xml", &xmlFile, "xml file with solver parameters for a 2x2 blocked NS example");

  switch (clp.parse(argc,argv)) {
  case CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case CommandLineProcessor::PARSE_ERROR:
  case CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  RCP<FancyOStream> out = fancyOStream(rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Time myTime("global");
  TimeMonitor MM(myTime);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  // custom parameters
  LO maxLevels = 2;   // TODO: singular system if MaxLevels > 2?

  GO maxCoarseSize=1; //FIXME clp doesn't like long long int

  int globalNumDofs = 8898;  // used for the maps
  int nDofsPerNode = 3;	     // used for generating the fine level null-space

  // build strided maps
  // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  /////////////////////////////////////// build strided maps
  // build strided maps:
  // xstridedfullmap: full map (velocity and pressure dof gids), continous
  // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
  // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
  RCP<const StridedMap> xstridedfullmap = StridedMapFactory::Build(lib,globalNumDofs,0,stridingInfo,comm,-1);
  RCP<const StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap,0);
  RCP<const StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap,1);

  /////////////////////////////////////// transform Xpetra::Map objects to Epetra
  // this is needed for AztecOO
  const RCP<const Epetra_Map> fullmap = rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
  RCP<const Epetra_Map>       velmap  = rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
  RCP<const Epetra_Map>       premap  = rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

  /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

  // read in problem
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  *out << "Reading matrix market file" << std::endl;
  EpetraExt::MatrixMarketFileToCrsMatrix("A5932_re1000.txt",*fullmap,*fullmap,*fullmap,ptrA);
  EpetraExt::MatrixMarketFileToVector("b5932_re1000.txt",*fullmap,ptrf);
  //EpetraExt::MatrixMarketFileToCrsMatrix("/home/wiesner/promotion/trilinos/fc16-debug/packages/muelu/test/navierstokes/A5932_re1000.txt",*fullmap,*fullmap,*fullmap,ptrA);
  //EpetraExt::MatrixMarketFileToVector("/home/wiesner/promotion/trilinos/fc16-debug/packages/muelu/test/navierstokes/b5932_re1000.txt",*fullmap,ptrf);
  //EpetraExt::MatrixMarketFileToCrsMatrix("/home/tobias/promotion/trilinos/fc17-dyn/packages/muelu/test/navierstokes/A5932_re1000.txt",*fullmap,*fullmap,*fullmap,ptrA);
  //EpetraExt::MatrixMarketFileToVector("/home/tobias/promotion/trilinos/fc17-dyn/packages/muelu/test/navierstokes/b5932_re1000.txt",*fullmap,ptrf);
  RCP<Epetra_CrsMatrix> epA = rcp(ptrA);
  RCP<Epetra_Vector> epv = rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = rcp(ptrNS);

  /////////////////////////////////////// split system into 2x2 block system

  *out << "Split matrix into 2x2 block matrix" << std::endl;

  // split fullA into A11,..., A22
  RCP<Epetra_CrsMatrix> A11;
  RCP<Epetra_CrsMatrix> A12;
  RCP<Epetra_CrsMatrix> A21;
  RCP<Epetra_CrsMatrix> A22;

  if(SplitMatrix2x2(epA,*velmap,*premap,A11,A12,A21,A22)==false)
    *out << "Problem with splitting matrix"<< std::endl;

  /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

  // build Xpetra objects from Epetra_CrsMatrix objects
  RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA11 = rcp(new Xpetra::EpetraCrsMatrix(A11));
  RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA12 = rcp(new Xpetra::EpetraCrsMatrix(A12));
  RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA21 = rcp(new Xpetra::EpetraCrsMatrix(A21));
  RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA22 = rcp(new Xpetra::EpetraCrsMatrix(A22));

  /////////////////////////////////////// generate MapExtractor object

  std::vector<RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
  xmaps.push_back(xstridedvelmap);
  xmaps.push_back(xstridedpremap);

  RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xstridedfullmap,xmaps);

  /////////////////////////////////////// build blocked transfer operator
  // using the map extractor
  RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_extractor,10));
  bOp->setMatrix(0,0,xA11);
  bOp->setMatrix(0,1,xA12);
  bOp->setMatrix(1,0,xA21);
  bOp->setMatrix(1,1,xA22);

  bOp->fillComplete();

  //////////////////////////////////////// prepare setup
  ParameterListInterpreter mueLuFactory(xmlFile, *comm);


  RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
  H->setDefaultVerbLevel(VERB_HIGH);
  H->SetMaxCoarseSize(maxCoarseSize);

  RCP<MueLu::Level> Finest = H->GetLevel(0);
  Finest->setDefaultVerbLevel(VERB_HIGH);
  Finest->Set("A",           rcp_dynamic_cast<Matrix>(bOp));


  ////////////////////////////////////////// prepare null space for A11
  RCP<MultiVector> nullspace11 = MultiVectorFactory::Build(xstridedvelmap, 2);  // this is a 2D standard null space

  for (int i=0; i<nDofsPerNode-1; ++i) {
    ArrayRCP<Scalar> nsValues = nullspace11->getDataNonConst(i);
    int numBlocks = nsValues.size() / (nDofsPerNode - 1);
    for (int j=0; j< numBlocks; ++j) {
      nsValues[j*(nDofsPerNode - 1) + i] = 1.0;
    }
  }

  Finest->Set("Nullspace1",nullspace11);

  ////////////////////////////////////////// prepare null space for A22
  RCP<MultiVector> nullspace22 = MultiVectorFactory::Build(xstridedpremap, 1);  // this is a 2D standard null space
  ArrayRCP<Scalar> nsValues22 = nullspace22->getDataNonConst(0);
  for (int j=0; j< nsValues22.size(); ++j) {
    nsValues22[j] = 1.0;
  }

  Finest->Set("Nullspace2",nullspace22);

  /////////////////////////////////// BEGIN setup

  mueLuFactory.SetupHierarchy(*H);

  ///////////////////////////////////// END setup

  *out << std::endl;

  RCP<MultiVector> xLsg = MultiVectorFactory::Build(xstridedfullmap,1);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    // Epetra_Vector -> Xpetra::Vector
    RCP<Vector> xRhs = rcp(new Xpetra::EpetraVector(epv));

    // calculate initial (absolute) residual
    Array<ScalarTraits<SC>::magnitudeType> norms(1);
    xRhs->norm2(norms);
    *out << "||x_0|| = " << norms[0] << std::endl;

    // apply ten multigrid iterations
    H->Iterate(*xRhs,100,*xLsg);

    // calculate and print residual
    RCP<MultiVector> xTmp = MultiVectorFactory::Build(xstridedfullmap,1);
    bOp->apply(*xLsg,*xTmp,NO_TRANS,(SC)1.0,(SC)0.0);
    xRhs->update((SC)-1.0,*xTmp,(SC)1.0);
    xRhs->norm2(norms);
    *out << "||r|| = " << norms[0] << std::endl;

  }

  // TODO: don't forget to add Aztec as prerequisite in CMakeLists.txt!
  //
  // Solve Ax = b using AMG as a preconditioner in AztecOO
  //
  {
    RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
    X->PutScalar(0.0);
    Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

    AztecOO aztecSolver(epetraProblem);
    aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

    MueLu::EpetraOperator aztecPrec(H);
    aztecSolver.SetPrecOperator(&aztecPrec);

    int maxIts = 50;
    double tol = 1e-8;

    aztecSolver.Iterate(maxIts, tol);
  }

   return EXIT_SUCCESS;
}
