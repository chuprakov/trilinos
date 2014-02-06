/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/


#ifndef CG_SOLVE_FILE_HPP
#define CG_SOLVE_FILE_HPP
#include "Tpetra_ConfigDefs.hpp"
#include "Kokkos_ConfigDefs.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLPerfTestArchive.hpp>
#include <Teuchos_Array.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include <TpetraUtils_MatrixGenerator.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <algorithm>
#include <functional>


#include <impl/Kokkos_Timer.hpp>

struct result_struct {
  double addtime,dottime,matvectime,final_residual;
  int niters;
  result_struct(double add, double dot, double matvec,int niter,double res):
    addtime(add),dottime(dot),matvectime(matvec),
    niters(niter),final_residual(res) {};
};


template<class Node>
Teuchos::XMLTestNode machine_configuration(Node node);

template<class Node>
Teuchos::XMLTestNode machine_configuration(Node node) {
  Teuchos::XMLTestNode config = Teuchos::PerfTest_MachineConfig();
  config.addString("KokkosNodeType",node->name());
  return config;
}


template<class CrsMatrix>
Teuchos::XMLTestNode test_entry(
    const std::string& filename_matrix,
    const std::string& filename_vector,
    int nsize, int mpi_ranks, int teams, int threads,
    Teuchos::RCP<CrsMatrix> A,
    result_struct results,int max_iters,int res_tolerance,
    double tol_small, double tol_large) {
  Teuchos::XMLTestNode configuration("TestConfiguration");
  configuration.addInt("MPI_Ranks",mpi_ranks);
  configuration.addInt("Teams",teams);
  configuration.addInt("Threads",threads);
  if(filename_matrix.empty()) {
    std::ostringstream strs;
    strs << "MiniFE-Generated " << nsize;
    configuration.addString("Matrix_File",strs.str());
  } else
    configuration.addString("Matrix_File",filename_matrix);

  if(filename_vector.empty()) {
    std::ostringstream strs;
    strs << "MiniFE-Generated " << nsize;
    configuration.addString("Vector_File",strs.str());
  } else
    configuration.addString("Vector_File",filename_vector);

  configuration.addInt("Matrix_Rows",A->getGlobalNumRows());
  configuration.addInt("Matrix_Cols",A->getGlobalNumCols());
  configuration.addInt("Matrix_NNZ",A->getGlobalNumEntries());
  configuration.addInt("Max_Iterations",max_iters);
  configuration.addInt("Tolerance",res_tolerance);


  Teuchos::XMLTestNode times("TestResults");
  times.addValueTolerance("Time_Add", Teuchos::ValueTolerance(results.addtime,tol_large));
  times.addValueTolerance("Time_Dot", Teuchos::ValueTolerance(results.dottime,tol_large));
  times.addValueTolerance("Time_MatVec", Teuchos::ValueTolerance(results.matvectime,tol_large));
  times.addValueTolerance("Time_CGSolve", Teuchos::ValueTolerance(results.matvectime+results.addtime+results.dottime,tol_small));
  times.addValueTolerance("Result_Iterations",Teuchos::ValueTolerance(results.niters,
                                                             results.niters>0?results.niters-1:0,
                                                             results.niters+1));
  times.addValueTolerance("Result_Final_Residual",Teuchos::ValueTolerance(results.final_residual,tol_small));

  //times.addString("Result_Residual", ValueTolerance(atof(argc[8]),4,6).as_string());

  Teuchos::XMLTestNode test("TPetra::CG-Solve");
  Teuchos::XMLTestNode entry("TestEntry");
  entry.addChild(configuration);
  entry.addChild(times);
  test.addChild(entry);
  return test;
}


template<class CrsMatrix, class Vector>
result_struct cg_solve(Teuchos::RCP<CrsMatrix> A, Teuchos::RCP<Vector> b, Teuchos::RCP<Vector> x, int myproc) {
  typedef double ScalarType;
  typedef double magnitude_type;
  typedef typename CrsMatrix::local_ordinal_type LocalOrdinalType;
  Teuchos::RCP<Vector> r,p,Ap;
  int max_iter=200;
  double tolerance = 1e-8;
  r = Tpetra::createVector<ScalarType>(A->getRangeMap());
  p = Tpetra::createVector<ScalarType>(A->getRangeMap());
  Ap = Tpetra::createVector<ScalarType>(A->getRangeMap());

  int length = r->getLocalLength();
  for(int i = 0;i<length;i++) {
    x->replaceLocalValue(i,0);
    r->replaceLocalValue(i,1);
    Ap->replaceLocalValue(i,1);
  }

  magnitude_type normr = 0;
  magnitude_type rtrans = 0;
  magnitude_type oldrtrans = 0;

  LocalOrdinalType print_freq = max_iter/10;
  if (print_freq>50) print_freq = 50;
  if (print_freq<1)  print_freq = 1;

  double dottime = 0;
  double addtime = 0;
  double matvectime = 0;

  Kokkos::Impl::Timer timer;
  p->update(1.0,*x,0.0,*x,0.0);
  addtime += timer.seconds(); timer.reset();


  A->apply(*p, *Ap);
  matvectime += timer.seconds(); timer.reset();

  r->update(1.0,*b,-1.0,*Ap,0.0);
  addtime += timer.seconds(); timer.reset();

  rtrans = r->dot(*r);
  dottime += timer.seconds(); timer.reset();

  normr = std::sqrt(rtrans);

  if (myproc == 0) {
    std::cout << "Initial Residual = "<< normr << std::endl;
  }

  magnitude_type brkdown_tol = std::numeric_limits<magnitude_type>::epsilon();
  LocalOrdinalType k;
  for(k=1; k <= max_iter && normr > tolerance; ++k) {
    if (k == 1) {
      p->update(1.0,*r,0.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    else {
      oldrtrans = rtrans;
      rtrans = r->dot(*r);
      dottime += timer.seconds(); timer.reset();
      magnitude_type beta = rtrans/oldrtrans;
      p->update(beta,*p,1.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    normr = std::sqrt(rtrans);
    if (myproc == 0 && (k%print_freq==0 || k==max_iter)) {
      std::cout << "Iteration = "<<k<<"   Residual = "<<normr<<std::endl;
    }

    magnitude_type alpha = 0;
    magnitude_type p_ap_dot = 0;
    A->apply(*p, *Ap);
    matvectime += timer.seconds(); timer.reset();
    p_ap_dot = Ap->dot(*p);
    dottime += timer.seconds(); timer.reset();

   if (p_ap_dot < brkdown_tol) {
      if (p_ap_dot < 0 ) {
        std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"<<std::endl;
        return result_struct(0,0,0,0,0);
      }
      else brkdown_tol = 0.1 * p_ap_dot;
    }
    alpha = rtrans/p_ap_dot;


    x->update(1.0,*x,alpha,*p,0.0);
    r->update(1.0,*r,-alpha,*Ap,0.0);
    addtime += timer.seconds(); timer.reset();

  }
  rtrans = r->dot(*r);

  normr = std::sqrt(rtrans);


  return result_struct(addtime,dottime,matvectime,k-1,normr);
}

template<class Node>
int run(int argc, char *argv[]) {
  typedef double                                                  Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType            Magnitude;
  typedef int                                                     Ordinal;

  typedef Tpetra::MpiPlatform<Node>                            Platform;
  typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal,Node>          CrsMatrix;
  using Teuchos::RCP;
  using Teuchos::tuple;


  //
  // Get example parameters from command-line processor
  //
  bool printMatrix = false;
  bool verbose = false;
  int niters = 100;
  int numthreads = 1;
  int numteams = 1;
  int nsize = 20;
  Magnitude tolerance = 1.0e-2;
  std::string filename;
  std::string filename_vector;
  std::string testarchive("Tpetra_PerformanceTests.xml");
  std::string hostname;

  double tol_small = 0.05;
  double tol_large = 0.10;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("numthreads",&numthreads,"Number of threads per thread team.");
  cmdp.setOption("numteams",&numteams,"Number of thread teams.");
  cmdp.setOption("hostname",&hostname,"Override of hostname for PerfTest entry.");
  cmdp.setOption("testarchive",&testarchive,"Set filename for Performance Test archive.");
  cmdp.setOption("filename",&filename,"Filename for test matrix.");
  cmdp.setOption("filename_vector",&filename_vector,"Filename for test matrix vector.");
  cmdp.setOption("tolerance",&tolerance,"Relative residual tolerance used for solver.");
  cmdp.setOption("iterations",&niters,"Maximum number of iterations.");
  cmdp.setOption("printMatrix","noPrintMatrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("size",&nsize,"Generate miniFE matrix with X^3 elements.");
  cmdp.setOption("tol_small",&tol_small,"Tolerance for total CG-Time and final residual.");
  cmdp.setOption("tol_large",&tol_small,"Tolerance for individual times.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return int(-1);
  }

  int myRank = 0;
#ifdef HAVE_MPI
  (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
#endif // HAVE_MPI

  int device = myRank%2;
  if (device>0) device++;
  int verboseint = verbose?1:0;
  Teuchos::ParameterList params;
  params.set("Num Threads",numthreads,"Number of Threads per Threadteam");
  params.set("Num Teams",numteams,"Number of Threadteams");
  params.set("Verbose",verboseint,"Verbose output");
  params.set("Device",device,"Device Number");



  //
  // Get the communicator and node
  //
  Node anode(params);
  RCP<Node>  node(&anode,false);

  Platform platform(node);
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();

  /*Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();*/
  //const int myRank = comm->getRank();
  //
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl << std::endl;
    std::cout << "Comm info: " << *comm;
  }


  // Read Tpetra::CrsMatrix from file
  //
  RCP<CrsMatrix> A;
  if(!filename.empty())
    A = Tpetra::MatrixMarket::Reader<CrsMatrix>::readSparseFile(filename,comm,node);
  else
    A = Tpetra::Utils::MatrixGenerator<CrsMatrix>::generate_miniFE_matrix(nsize,comm,node);
  if (printMatrix) {
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if (verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node> Vector;

  if ( A->getRangeMap() != A->getDomainMap() ) {
    throw std::runtime_error("TpetraExamples::powerMethod(): operator must have domain and range maps that are equivalent.");
  }
  // create three vectors, fill z with random numbers
  Teuchos::RCP<Vector> b, x;
  RCP<const typename CrsMatrix::map_type> map = A->getRangeMap();

  if(nsize<0)
    b = Tpetra::MatrixMarket::Reader<CrsMatrix>::readVectorFile(filename_vector,comm,node,map);
  else
    b = Tpetra::Utils::MatrixGenerator<CrsMatrix>::generate_miniFE_vector (nsize,comm,node);

  x = Tpetra::createVector<Scalar>(A->getRangeMap());

  result_struct results = cg_solve(A,b,x,myRank);
  if (myRank == 0) {
    Teuchos::XMLTestNode machine_config = machine_configuration(node);
    Teuchos::XMLTestNode test = test_entry(filename,filename_vector,nsize,
                                           comm->getSize(),numteams,numthreads,
                                           A,results,niters,tolerance,tol_small,tol_large);
    Teuchos::PerfTestResult comparison_result=Teuchos::PerfTest_CheckOrAdd_Test(machine_config,test,testarchive,hostname);
    switch (comparison_result) {
      case  Teuchos::PerfTestPassed: std::cout << "PASSED" << std::endl; break;
      case  Teuchos::PerfTestFailed: std::cout << "FAILED" << std::endl; break;
      case  Teuchos::PerfTestNewMachine: std::cout << "PASSED. Adding new machine entry." << std::endl; break;
      case  Teuchos::PerfTestNewConfiguration: std::cout << "PASSED. Adding new machine configuration." << std::endl; break;
      case  Teuchos::PerfTestNewTest: std::cout << "PASSED. Adding new test entry." << std::endl; break;
      case  Teuchos::PerfTestNewTestConfiguration: std::cout << "PASSED. Adding new test entry configuration." << std::endl; break;
      case  Teuchos::PerfTestUpdatedTest: std::cout << "PASSED. Updating test entry." << std::endl; break;
    }
    if(verbose) std::cout << test << std::endl;
  }

  return 0;
}

#endif
