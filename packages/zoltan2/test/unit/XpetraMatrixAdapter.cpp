// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of Zoltan2::XpetraCrsMatrixInput 

#include <string>


#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <UserInputForTests.hpp>

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

typedef double scalar_t;
typedef int lno_t;
typedef int gno_t;
typedef Zoltan2::default_node_t node_t;

typedef UserInputForTests<scalar_t, lno_t, gno_t> uinput_t;
typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tmatrix_t;
typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;
typedef Epetra_CrsMatrix ematrix_t;

template <typename L, typename G>
  void printMatrix(RCP<const Comm<int> > &comm, L nrows,
    const G *rowIds, const L *offsets, const G *colIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (L i=0; i < nrows; i++){
        std::cout << " row " << rowIds[i] << ": ";
        for (L j=offsets[i]; j < offsets[i+1]; j++){
          std::cout << colIds[j] << " ";
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
  Zoltan2::XpetraCrsMatrixInput<User> &ia, tmatrix_t &M)
{
  typedef typename Zoltan2::InputTraits<User>::scalar_t S;
  typedef typename Zoltan2::InputTraits<User>::lno_t L;
  typedef typename Zoltan2::InputTraits<User>::gno_t G;

  RCP<const Comm<int> > comm = M.getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getLocalNumRows() != M.getNodeNumRows())
    fail = 4;

  if (!fail && ia.getGlobalNumRows() != M.getGlobalNumRows())
    fail = 5;

  if (M.getNodeNumRows()){
    if (!fail && ia.getLocalNumColumns() != M.getNodeNumCols())
      fail = 6;
  }

  if (!fail && ia.getGlobalNumColumns() != M.getGlobalNumCols())
    fail = 7;

  gfail = globalFail(comm, fail);

  const G *rowIds=NULL, *colIds=NULL;
  const L *offsets=NULL;
  size_t nrows=0;

  if (!gfail){

    nrows = ia.getRowListView(rowIds, offsets, colIds);

    if (nrows != M.getNodeNumRows())
      fail = 8;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printMatrix<L, G>(comm, nrows, rowIds, offsets, colIds);
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

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra matrices for testing.

  RCP<uinput_t> uinput;

  try{
    uinput = 
      rcp(new uinput_t(std::string("../data/simple.mtx"), comm));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tmatrix_t> tM;     // original matrix (for checking)
  RCP<tmatrix_t> newM;   // migrated matrix

  tM = uinput->getTpetraCrsMatrix();
  size_t nrows = tM->getNodeNumRows();
  Teuchos::ArrayView<const gno_t> rowGids = 
    tM->getRowMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.

  typedef Zoltan2::IdentifierMap<gno_t, lno_t, gno_t> idmap_t;
  typedef Zoltan2::PartitioningSolution<scalar_t, gno_t, lno_t, gno_t> soln_t;

  RCP<const Zoltan2::Environment> env = Zoltan2::getDefaultEnvironment();

  ArrayRCP<const gno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, gidArray));

  int weightDim = 1;

  scalar_t *imbal = new scalar_t [weightDim];
  imbal[0] = 1.0;
  ArrayRCP<scalar_t> metric(imbal, 0, 1, true);

  size_t *p = new size_t [nvtx];
  memset(p, 0, sizeof(size_t) * nvtx);
  ArrayRCP<size_t> solnParts(p, 0, nvtx, true);

  soln_t solution(env, idMap, weightDim);

  solution.setParts(rowGids, solnParts, metric);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsMatrix
  if (!gfail){ 
    RCP<const tmatrix_t> ctM = rcp_const_cast<const tmatrix_t>(tM);
    RCP<Zoltan2::XpetraCrsMatrixInput<tmatrix_t> > tMInput;
  
    try {
      tMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixInput<tmatrix_t>(ctM));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraCrsMatrixInput ")+e.what(), 1);
    }
  
    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsMatrix" << std::endl;
    
    fail = verifyInputAdapter<tmatrix_t>(*tMInput, *tM);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      tmatrix_t *mMigrate = NULL;
      try{
        tMInput->applyPartitioningSolution(*tM, mMigrate, solution);
        newM = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const tmatrix_t> cnewM = rcp_const_cast<const tmatrix_t>(newM);
        RCP<Zoltan2::XpetraCrsMatrixInput<tmatrix_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsMatrixInput<tmatrix_t>(cnewM));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraCrsMatrixInput 2 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Tpetra::CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<tmatrix_t>(*newInput, *newM);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::CrsMatrix
  if (!gfail){ 
    RCP<xmatrix_t> xM = uinput->getXpetraCrsMatrix();
    RCP<const xmatrix_t> cxM = rcp_const_cast<const xmatrix_t>(xM);
    RCP<Zoltan2::XpetraCrsMatrixInput<xmatrix_t> > xMInput;
  
    try {
      xMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixInput<xmatrix_t>(cxM));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraCrsMatrixInput 3 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Xpetra::CrsMatrix" << std::endl;
    }
    fail = verifyInputAdapter<xmatrix_t>(*xMInput, *tM);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      xmatrix_t *mMigrate =NULL;
      try{
        xMInput->applyPartitioningSolution(*xM, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const xmatrix_t> cnewM(mMigrate);
        RCP<Zoltan2::XpetraCrsMatrixInput<xmatrix_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraCrsMatrixInput<xmatrix_t>(cnewM));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraCrsMatrixInput 4 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Xpetra::CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<xmatrix_t>(*newInput, *newM);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Epetra_CrsMatrix
  if (!gfail){ 
    RCP<ematrix_t> eM = uinput->getEpetraCrsMatrix();
    RCP<const ematrix_t> ceM = rcp_const_cast<const ematrix_t>(eM);
    RCP<Zoltan2::XpetraCrsMatrixInput<ematrix_t> > eMInput;
  
    try {
      eMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixInput<ematrix_t>(ceM));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraCrsMatrixInput 5 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Epetra_CrsMatrix" << std::endl;
    }
    fail = verifyInputAdapter<ematrix_t>(*eMInput, *tM);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      ematrix_t *mMigrate =NULL;
      try{
        eMInput->applyPartitioningSolution(*eM, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const ematrix_t> cnewM(mMigrate, true);
        RCP<Zoltan2::XpetraCrsMatrixInput<ematrix_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraCrsMatrixInput<ematrix_t>(cnewM));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraCrsMatrixInput 6 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Epetra_CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<ematrix_t>(*newInput, *newM);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
