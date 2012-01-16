// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Test for Zoltan2::BasicVectorInput 

#include <Zoltan2_BasicVectorInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;

typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> userTypes_t;

int checkBasicVector(
  Zoltan2::BasicVectorInput<userTypes_t> *ia, int len, int glen,
  gno_t *ids, int mvdim, scalar_t **values, int *valueStrides,
  int wdim, scalar_t **weights, int *weightStrides)
{
  int fail = 0;
  bool strideOne = false;

  if (valueStrides == NULL) strideOne = true;

  if (ia->getNumberOfVectors() != mvdim)
    fail = 100;

  if (!fail && ia->getNumberOfWeights() != wdim)
    fail = 101;

  if (!fail && ia->getLocalLength() != len)
    fail = 102;

  if (!fail && ia->getGlobalLength() != glen)
    fail = 103;

  for (int v=0; !fail && v < mvdim; v++){
    const gno_t *idList;
    const scalar_t *vals;
    int correctStride = (strideOne ? 1 : valueStrides[v]);
    int stride;

    size_t nvals = ia->getVector(v, idList, vals, stride);

    if (nvals != len*stride)
      fail = 104;

    if (!fail && stride != correctStride)
      fail = 105;

    for (int i=0; !fail && i < len; i++){
// TODO fix values check
//      if (vals[stride*i] != values[v][correctStride*i])
//        fail = 106;

      if (!fail && idList[i] != ids[i])
        fail = 107;
    }
  }

  for (int w=0; !fail && w < wdim; w++){
    const scalar_t *wgts;
    int stride;

    size_t nvals = ia->getVectorWeights(w, wgts, stride);

    if (nvals != len*stride)
      fail = 108;

    if (!fail && stride != weightStrides[w])
      fail = 109;

    for (int i=0; !fail && i < len; i++){
      if (wgts[stride*i] != weights[w][weightStrides[w]*i])
        fail = 110;
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
  int fail = 0;

  // Create a single vector and a strided multi-vector with 
  // strided multi-weights.

  lno_t numLocalIds = 10;
  gno_t *myIds = new gno_t [numLocalIds];
  gno_t base = rank * numLocalIds;
  
  int wdim = 2;
  scalar_t *weights = new scalar_t [numLocalIds*wdim];
  int *weightStrides = new int [wdim];
  scalar_t **weightPtrs = new scalar_t * [wdim];

  int mvdim = 3;
  scalar_t *v_values= new scalar_t [numLocalIds];
  scalar_t *mv_values= new scalar_t [mvdim * numLocalIds];
  int *valueStrides = new int [mvdim];
  scalar_t **valuePtrs = new scalar_t * [mvdim];

  for (lno_t i=0; i < numLocalIds; i++){
    myIds[i] = base+i;

    for (int w=0; w < wdim; w++)
      weights[w*numLocalIds + i] = w + 1 + nprocs - rank;

    v_values[i] = numLocalIds-i;

    for (int v=0; v < mvdim; v++)
      mv_values[i*mvdim + v] = (v+1) * (nprocs-rank) / (i+1);
  }

  for (int w=0; w < wdim; w++){
    weightStrides[w] = 1;
    weightPtrs[w] = weights + numLocalIds*w;
  }

  for (int v=0; v < mvdim; v++){
    valueStrides[v] = mvdim;
    valuePtrs[v] = mv_values + v;
  }

  RCP<Zoltan2::BasicVectorInput<userTypes_t> > ia;

  // A Zoltan2::BasicVectorInput object with one vector and no weights

  try{
   ia = rcp(new Zoltan2::BasicVectorInput<userTypes_t>(numLocalIds, myIds,
     v_values, 1, 0, NULL, NULL));
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 1", fail);

  fail = checkBasicVector(ia.getRawPtr(), numLocalIds, numLocalIds*nprocs,
    myIds, 1, valuePtrs, NULL, 0, NULL, NULL);

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 1", fail);

  // A Zoltan2::BasicVectorInput object with one vector and weights

  try{
   ia = rcp(new Zoltan2::BasicVectorInput<userTypes_t>(numLocalIds, myIds,
     v_values, 1, wdim, weightPtrs, weightStrides));
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 2", fail);

  fail = checkBasicVector(ia.getRawPtr(), numLocalIds, numLocalIds*nprocs,
    myIds, 1, valuePtrs, NULL, wdim, weightPtrs, weightStrides);

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 2", fail);

  // A Zoltan2::BasicVectorInput object with a multivector and no weights

  try{
   ia = rcp(new Zoltan2::BasicVectorInput<userTypes_t>(mvdim,
    numLocalIds, myIds, valuePtrs, valueStrides,
    0, NULL, NULL));
    
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 3", fail);

  fail = checkBasicVector(ia.getRawPtr(), numLocalIds, numLocalIds*nprocs,
    myIds, mvdim, valuePtrs, valueStrides, 0, NULL, NULL);

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 3", fail);

  // A Zoltan2::BasicVectorInput object with a multivector with weights

  try{
   ia = rcp(new Zoltan2::BasicVectorInput<userTypes_t>(mvdim,
    numLocalIds, myIds, valuePtrs, valueStrides,
    wdim, weightPtrs, weightStrides));
    
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);

  fail = checkBasicVector(ia.getRawPtr(), numLocalIds, numLocalIds*nprocs,
    myIds, mvdim, valuePtrs, valueStrides,
    wdim, weightPtrs, weightStrides);

  TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 4", fail);

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return fail;
}

