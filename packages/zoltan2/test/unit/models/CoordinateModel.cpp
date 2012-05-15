// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Testing of CoordinateModel
//

#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <set>

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_CrsMatrix.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;

void testCoordinateModel(std::string &fname, int weightDim,
  const RCP<const Comm<int> > &comm, bool consecutiveIds)
{
  int fail = 0, gfail = 0;

  //////////////////////////////////////////////////////////////
  // Input data
  //////////////////////////////////////////////////////////////

  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> mv_t;

  RCP<UserInputForTests> uinput;

  try{
    uinput = rcp(new UserInputForTests(testDataFilePath, fname, comm, true));
  }
  catch(std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "input constructor", 1);

  RCP<mv_t> coords;

  try{
    coords = uinput->getCoordinates();
  }
  catch(std::exception &e){
    fail=2;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "getting coordinates", 1);

  int coordDim = coords->getNumVectors();

  TEST_FAIL_AND_EXIT(*comm, coordDim <= 3, "dim 3 at most", 1);

  const scalar_t *x=NULL, *y=NULL, *z=NULL;

  x = coords->getData(0).getRawPtr();
  if (coordDim > 1){
    y = coords->getData(1).getRawPtr();
    if (coordDim > 2)
      z = coords->getData(2).getRawPtr();
  }

  // Are these coordinates correct

  int nLocalIds = coords->getLocalLength();
  int nGlobalIds = coords->getGlobalLength();
  ArrayView<const gno_t> idList = coords->getMap()->getNodeElementList();

  Array<ArrayRCP<const scalar_t> > coordWeights(weightDim);

  for (int wdim=0; wdim < weightDim; wdim++){
    scalar_t *w = new scalar_t [nLocalIds];
    for (int i=0; i < nLocalIds; i++){
      w[i] = ((i%2) + 1) + wdim;
    }
    coordWeights[wdim] = Teuchos::arcp(w, 0, nLocalIds);
  }

  //////////////////////////////////////////////////////////////
  // Create a BasicCoordinateInput adapter object.
  //////////////////////////////////////////////////////////////

  typedef Zoltan2::BasicCoordinateInput<mv_t> ia_t;
  typedef Zoltan2::CoordinateInput<mv_t>      base_ia_t;

  RCP<ia_t> ia;

  if (weightDim == 0){   // use the simpler constructor
    try{
      ia = rcp(new ia_t(nLocalIds, idList.getRawPtr(), x, y, z));
    }
    catch(std::exception &e){
      fail=3;
    }
  }
  else{
    std::vector<const scalar_t *> values, weights;
    std::vector<int> valueStrides, weightStrides;  // default is 1
    values.push_back(x);
    if (y) {
      values.push_back(y);
      if (z) 
        values.push_back(z);
    }
    for (int wdim=0; wdim < weightDim; wdim++){
      weights.push_back(coordWeights[wdim].getRawPtr());
    }

    try{
      ia = rcp(new ia_t(nLocalIds, idList.getRawPtr(),
               values, valueStrides, weights, weightStrides));
    }
    catch(std::exception &e){
      fail=4;
    }
  }

  RCP<base_ia_t> base_ia = Teuchos::rcp_implicit_cast<base_ia_t>(ia);

  TEST_FAIL_AND_EXIT(*comm, !fail, "making input adapter", 1);

  //////////////////////////////////////////////////////////////
  // Create an CoordinateModel with this input
  //////////////////////////////////////////////////////////////

  typedef Zoltan2::StridedData<lno_t, scalar_t> input_t;
  typedef std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags_t;
  typedef Zoltan2::CoordinateModel<base_ia_t> model_t;
  modelFlags_t modelFlags;

  if (consecutiveIds)
    modelFlags.set(Zoltan2::IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);
  RCP<model_t> model;
  

  try{
    model = rcp(new model_t(base_ia.getRawPtr(), env, comm, modelFlags));
  }
  catch (std::exception &e){
    fail = 5;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "making model", 1);

  // Test the CoordinateModel interface

  if (model->getCoordinateDim() != coordDim)
    fail = 6;

  if (!fail && model->getLocalNumCoordinates() != size_t(nLocalIds))
    fail = 7;

  if (!fail && model->getGlobalNumCoordinates() != size_t(nGlobalIds))
    fail = 8;

  if (!fail && model->getCoordinateWeightDim() !=  weightDim)
    fail = 9;

  gfail = globalFail(comm, fail);

  if (gfail)
    printFailureCode(comm, fail);
  
  ArrayView<const gno_t> gids;
  ArrayView<input_t> xyz;
  ArrayView<input_t> wgts;
  
  model->getCoordinates(gids, xyz, wgts);

  if (!fail && gids.size() != nLocalIds)
    fail = 10;

  for (int i=0; !fail && i < nLocalIds; i++){
    if (gids[i] != idList[i])
      fail = 11;
  }

  if (!fail && wgts.size() != weightDim)
    fail = 12;

  const scalar_t *vals[3] = {x, y, z};

  for (int dim=0; !fail && dim < coordDim; dim++){
    for (int i=0; !fail && i < nLocalIds; i++){
      if (xyz[dim][i] != vals[dim][i])
        fail = 13;
    }
  }

  for (int wdim=0; !fail && wdim < weightDim; wdim++){
    for (int i=0; !fail && i < nLocalIds; i++){
      if (wgts[wdim][i] != coordWeights[wdim][i])
        fail = 14;
    }
  }

  if (!fail && consecutiveIds){
    bool inARow = Zoltan2::IdentifierTraits<gno_t>::areConsecutive(
      gids.getRawPtr(), nLocalIds);

    if (!inARow)
      fail = 15;
  }

  gfail = globalFail(comm, fail);

  if (gfail)
    printFailureCode(comm, fail);
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  string fname("commanche_dual");   // reader will seek coord file
  bool wishConsecutiveIds = true;

  if (rank == 0){
    std::cout << std::endl << fname;
    std::cout << ", without consecutive IDs, no weights" << std::endl;
  }
  testCoordinateModel(fname, 0, comm, !wishConsecutiveIds);

  if (rank == 0){
    std::cout << std::endl << fname;
    std::cout << ", with consecutive IDs, no weights" << std::endl;
  }
  testCoordinateModel(fname, 0, comm,  wishConsecutiveIds);

  if (rank == 0){
    std::cout << std::endl << fname;
    std::cout << ", without consecutive IDs, dim 1 weights" << std::endl;
  }
  testCoordinateModel(fname, 1, comm, !wishConsecutiveIds);

  if (rank == 0){
    std::cout << std::endl << fname;
    std::cout << ", with consecutive IDs, dim 2 weights " << std::endl;
  }
  testCoordinateModel(fname, 2, comm,  wishConsecutiveIds);

  if (rank==0) std::cout << "PASS" << std::endl;

  return 0;
}
