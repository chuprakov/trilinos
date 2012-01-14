// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of Zoltan2::XpetraCrsGraphInput

#include <string>

#include <Zoltan2_XpetraCrsGraphInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::Array;
using Teuchos::ArrayView;

typedef UserInputForTests<scalar_t, lno_t, gno_t> uinput_t;
typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tgraph_t;
typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
typedef Epetra_CrsGraph egraph_t;

template <typename L, typename G>
  void printGraph(RCP<const Comm<int> > &comm, L nvtx,
    const G *vtxIds, const L *offsets, const G *edgeIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (L i=0; i < nvtx; i++){
        std::cout << " vertex " << vtxIds[i] << ": ";
        for (L j=offsets[i]; j < offsets[i+1]; j++){
          std::cout << edgeIds[j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout.flush();
    }
    comm->barrier();
  }
  comm->barrier();
}

template <typename User>
int verifyInputAdapter(
  Zoltan2::XpetraCrsGraphInput<User> &ia, tgraph_t &graph)
{
  RCP<const Comm<int> > comm = graph.getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getLocalNumVertices() != graph.getNodeNumRows())
    fail = 4;

  if (!fail && ia.getGlobalNumVertices() != graph.getGlobalNumRows())
    fail = 5;

  if (!fail && ia.getLocalNumEdges() != graph.getNodeNumEntries())
      fail = 6;

  if (!fail && ia.getGlobalNumEdges() != graph.getGlobalNumEntries())
    fail = 7;

  gfail = globalFail(comm, fail);

  const gno_t *vtxIds=NULL, *edgeIds=NULL;
  const lno_t *offsets=NULL;
  size_t nvtx=0;

  if (!gfail){

    nvtx = ia.getVertexListView(vtxIds, offsets, edgeIds);

    if (nvtx != graph.getNodeNumRows())
      fail = 8;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printGraph<lno_t, gno_t>(comm, nvtx, vtxIds, offsets, edgeIds);
    }
    else{
      if (!fail) fail = 10;
    }
  }
  return fail;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail=0;

  // Create an object that can give us test Tpetra, Xpetra
  // and Epetra graphs for testing.

  RCP<uinput_t> uinput;

  try{
    uinput =
      rcp(new uinput_t(std::string("../data/simple.mtx"), comm));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tgraph_t> tG;     // original graph (for checking)
  RCP<tgraph_t> newG;   // migrated graph

  tG = uinput->getTpetraCrsGraph();
  size_t nvtx = tG->getNodeNumRows();
  ArrayView<const gno_t> rowGids = tG->getRowMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.
  // Our solution just assigns all objects to part zero.

  typedef Zoltan2::IdentifierMap<tgraph_t> idmap_t;
  typedef Zoltan2::PartitioningSolution<tgraph_t> soln_t;

  RCP<const Zoltan2::Environment> env = Zoltan2::getDefaultEnvironment();

  ArrayRCP<const gno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, comm, gidArray));

  int weightDim = 1;

  float *imbal = new float [weightDim];
  imbal[0] = 1.0;
  ArrayRCP<float> metric(imbal, 0, 1, true);

  size_t *p = new size_t [nvtx];
  memset(p, 0, sizeof(size_t) * nvtx);
  ArrayRCP<size_t> solnParts(p, 0, nvtx, true);

  soln_t solution(env, comm, idMap, weightDim);

  solution.setParts(rowGids, solnParts, metric);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsGraph
  if (!gfail){
    RCP<const tgraph_t> ctG = rcp_const_cast<const tgraph_t>(tG);
    RCP<Zoltan2::XpetraCrsGraphInput<tgraph_t> > tGInput;

    try {
      tGInput =
        rcp(new Zoltan2::XpetraCrsGraphInput<tgraph_t>(ctG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphInput ")+e.what(), 1);
    }

    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsGraph" << std::endl;

    fail = verifyInputAdapter<tgraph_t>(*tGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      tgraph_t *mMigrate = NULL;
      try{
        tGInput->applyPartitioningSolution<tgraph_t>(*tG, mMigrate, solution);
        newG = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const tgraph_t> cnewG = rcp_const_cast<const tgraph_t>(newG);
        RCP<Zoltan2::XpetraCrsGraphInput<tgraph_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsGraphInput<tgraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphInput 2 ")+e.what(), 1);
        }

        if (rank==0){
          std::cout <<
           "Input adapter for Tpetra::CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<tgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::CrsGraph
  if (!gfail){
    RCP<xgraph_t> xG = uinput->getXpetraCrsGraph();
    RCP<const xgraph_t> cxG = rcp_const_cast<const xgraph_t>(xG);
    RCP<Zoltan2::XpetraCrsGraphInput<xgraph_t> > xGInput;

    try {
      xGInput =
        rcp(new Zoltan2::XpetraCrsGraphInput<xgraph_t>(cxG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphInput 3 ")+e.what(), 1);
    }

    if (rank==0){
      std::cout << "Input adapter for Xpetra::CrsGraph" << std::endl;
    }
    fail = verifyInputAdapter<xgraph_t>(*xGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      xgraph_t *mMigrate =NULL;
      try{
        xGInput->applyPartitioningSolution<tgraph_t>(*xG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const xgraph_t> cnewG(mMigrate);
        RCP<Zoltan2::XpetraCrsGraphInput<xgraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphInput<xgraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphInput 4 ")+e.what(), 1);
        }

        if (rank==0){
          std::cout <<
           "Input adapter for Xpetra::CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<xgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

#ifdef HAVE_EPETRA_DATA_TYPES
  /////////////////////////////////////////////////////////////
  // User object is Epetra_CrsGraph
  if (!gfail){
    RCP<egraph_t> eG = uinput->getEpetraCrsGraph();
    RCP<const egraph_t> ceG = rcp_const_cast<const egraph_t>(eG);
    RCP<Zoltan2::XpetraCrsGraphInput<egraph_t> > eGInput;

    try {
      eGInput =
        rcp(new Zoltan2::XpetraCrsGraphInput<egraph_t>(ceG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphInput 5 ")+e.what(), 1);
    }

    if (rank==0){
      std::cout << "Input adapter for Epetra_CrsGraph" << std::endl;
    }
    fail = verifyInputAdapter<egraph_t>(*eGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      egraph_t *mMigrate =NULL;
      try{
        eGInput->applyPartitioningSolution<tgraph_t>(*eG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const egraph_t> cnewG(mMigrate, true);
        RCP<Zoltan2::XpetraCrsGraphInput<egraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphInput<egraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphInput 6 ")+e.what(), 1);
        }

        if (rank==0){
          std::cout <<
           "Input adapter for Epetra_CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<egraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }
#endif

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
