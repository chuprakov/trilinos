// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef _ZOLTAN2_ALGSCOTCH_HPP_
#define _ZOLTAN2_ALGSCOTCH_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Scotch.hpp
//! \brief Parallel graph partitioning using Scotch.

namespace Zoltan2 {


#ifdef HAVE_ZOLTAN2_SCOTCH

// stdint.h for int64_t in scotch header
#include <stdint.h>
#ifndef HAVE_ZOLTAN2_MPI
#include "scotch.h"
#else
#include "ptscotch.h"
#endif

#ifdef SHOW_ZOLTAN2_SCOTCH_MEMORY
//
// Scotch keeps track of memory high water mark, but doesn't
// provide a way to get that number.  So add this function:
//   "size_t SCOTCH_getMemoryMax() { return memorymax;}"
// to src/libscotch/common_memory.c
//
// and this macro:
//   "#define HAVE_SCOTCH_GETMEMORYMAX
// to include/ptscotch.h
//
// and compile scotch with -DCOMMON_MEMORY_TRACE
//
#ifdef HAVE_SCOTCH_GETMEMORYMAX
  extern "C"{
    extern size_t SCOTCH_getMemoryMax();
  }
#else
#error "Either turn off ZOLTAN2_ENABLE_SCOTCH_MEMORY_REPORT in cmake configure, or see SHOW_ZOLTAN2_SCOTCH_MEMORY comment in Zoltan2_AlgScotch.hpp"
#endif  // HAVE_SCOTCH_GETMEMORYMAX
#endif  // SHOW_ZOLTAN2_SCOTCH_MEMORY



/*! Scotch partitioning method.
 *
 *  \param env  parameters for the problem and library configuration
 *  \param problemComm  the communicator for the problem
 *  \param model a graph
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).
 */



/////////////////////////////////////////////////////////////////////////////
//  Traits struct to handle conversions between gno_t/lno_t and SCOTCH_Num.
/////////////////////////////////////////////////////////////////////////////

// General case:  SCOTCH_Num and gno_t/lno_t (called zno_t here) differ.
template <typename zno_t>
struct SCOTCH_Num_Traits {

  static inline SCOTCH_Num ASSIGN_TO_SCOTCH_NUM(
    SCOTCH_Num &a,
    zno_t b,
    const RCP<const Environment> &env)
  {
    // Assign a = b; make sure SCOTCH_Num is large enough to accept zno_t.
    if (b <= SCOTCH_NUMMAX) a = b;
    else 
      env->localInputAssertion(__FILE__, __LINE__, 
       "Value too large for SCOTCH_Num, Rebuild Scotch with larger SCOTCH_Num",
       false, BASIC_ASSERTION);
    return a;
  }

  static inline void ASSIGN_SCOTCH_NUM_ARRAY(
    SCOTCH_Num **a,
    ArrayView<const zno_t> &b,
    const RCP<const Environment> &env)
  {
    // Allocate array a; copy b values into a.
    size_t size = b.size();
    if (size > 0) {
      *a = new SCOTCH_Num[size];
      for (size_t i = 0; i < size; i++) ASSIGN_TO_SCOTCH_NUM((*a)[i], b[i], env);
    }
    else {
      *a = NULL;
      // Note:  the Scotch manual says that if any rank has a non-NULL array,
      //        every process must have a non-NULL array.  In practice, 
      //        however, this condition is not needed for the arrays we use.
      //        For now, we'll set these arrays to NULL.  We could allocate
      //        a dummy value here if needed.  KDD 1/23/14
    }
  }

  static inline void DELETE_SCOTCH_NUM_ARRAY(SCOTCH_Num **a)
  {
    // Delete the copy made in ASSIGN_SCOTCH_NUM_ARRAY.
    delete [] *a;
  }
};


// Special case:  zno_t == SCOTCH_Num. No error checking or copies needed.
template <>
struct SCOTCH_Num_Traits<SCOTCH_Num> {
  static inline SCOTCH_Num ASSIGN_TO_SCOTCH_NUM(
    SCOTCH_Num &a,
    SCOTCH_Num b,
    const RCP<const Environment> &env)
  {
    a = b;
    return a;
  }
  static inline void ASSIGN_SCOTCH_NUM_ARRAY(
    SCOTCH_Num **a,
    ArrayView<const SCOTCH_Num> &b,
    const RCP<const Environment> &env)
  {
    if (b.size() > 0)
      *a = const_cast<SCOTCH_Num *> (b.getRawPtr());
    else {
      *a = NULL;
      // Note:  the Scotch manual says that if any rank has a non-NULL array,
      //        every process must have a non-NULL array.  In practice, 
      //        however, this condition is not needed for the arrays we use.
      //        For now, we'll set these arrays to NULL, because if we
      //        allocated a dummy value here, we'll have to track whether or
      //        not we can free it.  KDD 1/23/14
    }
  }
  static inline void DELETE_SCOTCH_NUM_ARRAY(SCOTCH_Num **a) { }
};

///////////////////////////////////////////////////////////////////////
// Now, the actual Scotch algorithm.
///////////////////////////////////////////////////////////////////////

template <typename Adapter>
class AlgPTScotch : public Algorithm<Adapter>
{
public:

  typedef GraphModel<typename Adapter::base_adapter_t> graphModel_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;

  /*! Scotch constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param model a graph
   *
   *  Preconditions: The parameters in the environment have been processed.
   *  TODO:  THIS IS A MINIMAL CONSTRUCTOR FOR NOW.
   *  TODO:  WHEN ADD SCOTCH ORDERING OR MAPPING, MOVE SCOTCH GRAPH CONSTRUCTION
   *  TODO:  TO THE CONSTRUCTOR SO THAT CODE MAY BE SHARED.
   */
  AlgPTScotch(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
#ifdef HAVE_ZOLTAN2_MPI
              MPI_Comm mpicomm__,
#endif
              const RCP<graphModel_t> &model__) :
    env(env__), problemComm(problemComm__), 
#ifdef HAVE_ZOLTAN2_MPI
    mpicomm(mpicomm__),
#endif
    model(model__)
  { }

  void partition(PartitioningSolution<Adapter> &);

private:

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpicomm;
#endif
  const RCP<GraphModel<typename Adapter::base_adapter_t> > model;

  void scale_weights(size_t n, StridedData<lno_t, scalar_t> &fwgts,
                     SCOTCH_Num *iwgts);
};


/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgPTScotch<Adapter>::partition(PartitioningSolution<Adapter> &solution)
{
  HELLO;

  size_t numGlobalParts = solution.getTargetGlobalNumberOfParts();

  SCOTCH_Num partnbr;
  SCOTCH_Num_Traits<size_t>::ASSIGN_TO_SCOTCH_NUM(partnbr, numGlobalParts, env);

#ifdef HAVE_ZOLTAN2_MPI
  int ierr = 0;
  int me = problemComm->getRank();

  const SCOTCH_Num  baseval = 0;  // Base value for array indexing.
                                  // GraphModel returns GNOs from base 0.

  SCOTCH_Strat stratstr;          // Strategy string
                                  // TODO:  Set from parameters
  SCOTCH_stratInit(&stratstr);

  // Allocate and initialize PTScotch Graph data structure.
  SCOTCH_Dgraph *gr = SCOTCH_dgraphAlloc();  // Scotch distributed graph
  ierr = SCOTCH_dgraphInit(gr, mpicomm);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphInit", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Get vertex info
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > xyz;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  size_t nVtx = model->getVertexList(vtxID, xyz, vwgts);
  SCOTCH_Num vertlocnbr;
  SCOTCH_Num_Traits<size_t>::ASSIGN_TO_SCOTCH_NUM(vertlocnbr, nVtx, env);
  SCOTCH_Num vertlocmax = vertlocnbr; // Assumes no holes in global nums.

  // Get edge info
  ArrayView<const gno_t> edgeIds;
  ArrayView<const int>   procIds;
  ArrayView<const lno_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;

  size_t nEdge = model->getEdgeList(edgeIds, procIds, offsets, ewgts);

  SCOTCH_Num edgelocnbr;
  SCOTCH_Num_Traits<size_t>::ASSIGN_TO_SCOTCH_NUM(edgelocnbr, nEdge, env);
  const SCOTCH_Num edgelocsize = edgelocnbr;  // Assumes adj array is compact.

  SCOTCH_Num *vertloctab;  // starting adj/vtx
  SCOTCH_Num_Traits<lno_t>::ASSIGN_SCOTCH_NUM_ARRAY(&vertloctab, offsets, env);

  SCOTCH_Num *edgeloctab;  // adjacencies
  SCOTCH_Num_Traits<gno_t>::ASSIGN_SCOTCH_NUM_ARRAY(&edgeloctab, edgeIds, env);

  // We don't use these arrays, but we need them as arguments to Scotch.
  SCOTCH_Num *vendloctab = NULL;  // Assume consecutive storage for adj
  SCOTCH_Num *vlblloctab = NULL;  // Vertex label array
  SCOTCH_Num *edgegsttab = NULL;  // Array for ghost vertices

  // Get weight info.
  SCOTCH_Num *velotab = NULL;  // Vertex weights
  SCOTCH_Num *edlotab = NULL;  // Edge weights

  int nVwgts = model->getNumWeightsPerVertex();
  int nEwgts = model->getNumWeightsPerEdge();
  if (nVwgts > 1 && me == 0) {
    std::cerr << "Warning:  NumWeightsPerVertex is " << nVwgts 
              << " but Scotch allows only one weight. "
              << " Zoltan2 will use only the first weight per vertex."
              << std::endl;
  }
  if (nEwgts > 1 && me == 0) {
    std::cerr << "Warning:  NumWeightsPerEdge is " << nEwgts 
              << " but Scotch allows only one weight. "
              << " Zoltan2 will use only the first weight per edge."
              << std::endl;
  }

  if (nVwgts) {
    velotab = new SCOTCH_Num[nVtx+1];  // +1 since Scotch wants all procs 
                                       // to have non-NULL arrays
    scale_weights(nVtx, vwgts[0], velotab);
  }

  if (nEwgts) {
    edlotab = new SCOTCH_Num[nEdge+1];  // +1 since Scotch wants all procs 
                                         // to have non-NULL arrays
    scale_weights(nEdge, ewgts[0], edlotab);
  }

  // Build PTScotch distributed data structure
  ierr = SCOTCH_dgraphBuild(gr, baseval, vertlocnbr, vertlocmax,
                            vertloctab, vendloctab, velotab, vlblloctab,
                            edgelocnbr, edgelocsize,
                            edgeloctab, edgegsttab, edlotab);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphBuild", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Create array for Scotch to return results in.
  ArrayRCP<part_t> partList(new part_t[nVtx], 0, nVtx,true);
  SCOTCH_Num *partloctab = NULL;
  if (nVtx && (sizeof(SCOTCH_Num) == sizeof(part_t))) {
    // Can write directly into the solution's memory
    partloctab = (SCOTCH_Num *) partList.getRawPtr();
  }
  else {
    // Can't use solution memory directly; will have to copy later.
    // Note:  Scotch does not like NULL arrays, so add 1 to always have non-null.
    //        ParMETIS has this same "feature."  See Zoltan bug 4299.
    partloctab = new SCOTCH_Num[nVtx+1];
  }

  // Get target part sizes
  float *partsizes = new float[numGlobalParts];
  if (!solution.criteriaHasUniformPartSizes(0))
    for (size_t i=0; i<numGlobalParts; i++)
      partsizes[i] = solution.getCriteriaPartSize(0, i);
  else
    for (size_t i=0; i<numGlobalParts; i++)
      partsizes[i] = 1.0 / float(numGlobalParts);

  // Allocate and initialize PTScotch target architecture data structure
  SCOTCH_Arch archdat;
  SCOTCH_archInit(&archdat);

  SCOTCH_Num velosum = 0;
  SCOTCH_dgraphSize (gr, &velosum, NULL, NULL, NULL);
  SCOTCH_Num *goalsizes = new SCOTCH_Num[partnbr];
  // TODO: The goalsizes are set as in Zoltan; not sure it is correct there 
  // or here.
  // It appears velosum is global NUMBER of vertices, not global total 
  // vertex weight.  I think we should use the latter.
  // Fix this when we add vertex weights.
  for (SCOTCH_Num i = 0; i < partnbr; i++)
    goalsizes[i] = ceil(velosum * partsizes[i]);
  delete [] partsizes;

  SCOTCH_archCmpltw(&archdat, partnbr, goalsizes);

  // Call partitioning; result returned in partloctab.
  ierr = SCOTCH_dgraphMap(gr, &archdat, &stratstr, partloctab);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphMap", 
    !ierr, BASIC_ASSERTION, problemComm);

  SCOTCH_archExit(&archdat);
  delete [] goalsizes;

  // TODO - metrics

#ifdef SHOW_ZOLTAN2_SCOTCH_MEMORY
  int me = env->comm_->getRank();
#endif

#ifdef HAVE_SCOTCH_ZOLTAN2_GETMEMORYMAX
  if (me == 0){
    size_t scotchBytes = SCOTCH_getMemoryMax();
    std::cout << "Rank " << me << ": Maximum bytes used by Scotch: ";
    std::cout << scotchBytes << std::endl;
  }
#endif

  // Clean up PTScotch
  SCOTCH_dgraphExit(gr);
  free(gr);
  SCOTCH_stratExit(&stratstr);

  // Load answer into the solution.

  if ((sizeof(SCOTCH_Num) != sizeof(part_t)) || (nVtx == 0)) {
    for (size_t i = 0; i < nVtx; i++) partList[i] = partloctab[i];
    delete [] partloctab;
  }

  ArrayRCP<const gno_t> gnos = arcpFromArrayView(vtxID);

  solution.setParts(gnos, partList, true);

  env->memory("Zoltan2-Scotch: After creating solution");

  // Clean up copies made due to differing data sizes.
  SCOTCH_Num_Traits<lno_t>::DELETE_SCOTCH_NUM_ARRAY(&vertloctab);
  SCOTCH_Num_Traits<gno_t>::DELETE_SCOTCH_NUM_ARRAY(&edgeloctab);

  if (nVwgts) delete [] velotab;
  if (nEwgts) delete [] edlotab;

#else // DO NOT HAVE_MPI

  // TODO:  Handle serial case with calls to Scotch.
  // TODO:  For now, assign everything to rank 0 and assume only one part.
  // TODO:  Can probably use the code above for loading solution,
  // TODO:  instead of duplicating it here.
  // TODO
  // TODO:  Actual logic should call Scotch when number of processes == 1.
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > xyz;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  size_t nVtx = model->getVertexList(vtxID, xyz, vwgts);

  ArrayRCP<part_t> partList(new part_t[nVtx], 0, nVtx, true);
  for (size_t i = 0; i < nVtx; i++) partList[i] = 0;

  ArrayRCP<const gno_t> gnos = arcpFromArrayView(vtxID);
  solution.setParts(gnos, partList, true);

#endif // DO NOT HAVE_MPI
}

/////////////////////////////////////////////////////////////////////////////
// Scale and round scalar_t weights (typically float or double) to 
// SCOTCH_Num (typically int or long).
// subject to sum(weights) <= max_wgt_sum.
// Only scale if deemed necessary.
//
// Note that we use ceil() instead of round() to avoid
// rounding to zero weights.
// Based on Zoltan's scale_round_weights, mode 1.

template <typename Adapter>
void AlgPTScotch<Adapter>::scale_weights(
  size_t n,
  StridedData<typename Adapter::lno_t, typename Adapter::scalar_t> &fwgts,
  SCOTCH_Num *iwgts
)
{
  const double INT_EPSILON = 1e-5;

  SCOTCH_Num nonint, nonint_local = 0;
  double sum_wgt, sum_wgt_local = 0.;
  double max_wgt, max_wgt_local = 0.;

  // Compute local sums of the weights 
  // Check whether all weights are integers
  for (size_t i = 0; i < n; i++) {
    double fw = double(fwgts[i]);
    if (!nonint_local){
      SCOTCH_Num tmp = (SCOTCH_Num) floor(fw + .5); /* Nearest int */
      if (fabs((double)tmp-fw) > INT_EPSILON) {
        nonint_local = 1;
      }
    }
    sum_wgt_local += fw;
    if (fw > max_wgt_local) max_wgt_local = fw;
  }

  Teuchos::reduceAll<int,int>(*problemComm, Teuchos::REDUCE_MAX, 1,
                              &nonint_local,  &nonint);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_SUM, 1,
                                 &sum_wgt_local, &sum_wgt);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_MAX, 1,
                                 &max_wgt_local, &max_wgt);

  double scale = 1.;
  const double max_wgt_sum = double(SCOTCH_NUMMAX/8);

  // Scaling needed if weights are not integers or weights' 
  // range is not sufficient
  if (nonint || (max_wgt <= INT_EPSILON) || (sum_wgt > max_wgt_sum)) {
    /* Calculate scale factor */
    if (sum_wgt != 0.) scale = max_wgt_sum/sum_wgt;
  }

  /* Convert weights to positive integers using the computed scale factor */
  for (size_t i = 0; i < n; i++)
    iwgts[i] = (SCOTCH_Num) ceil(double(fwgts[i])*scale);

}

////////////////////////////////////////////////////////////////////////
#else // DO NOT HAVE_ZOLTAN2_SCOTCH

// Error handling for when Scotch is requested
// but Zoltan2 not built with Scotch.

template <typename Adapter>
class AlgPTScotch : public Algorithm<Adapter>
{
public:
  AlgPTScotch(const RCP<const Environment> &env,
              const RCP<const Comm<int> > &problemComm,
#ifdef HAVE_ZOLTAN2_MPI
              MPI_Comm mpicomm,
#endif
              const RCP<GraphModel<typename Adapter::base_adapter_t> > &model
  )
  {
    throw std::runtime_error(
          "BUILD ERROR:  Scotch requested but not compiled into Zoltan2.\n"
          "Please set CMake flag Zoltan2_ENABLE_Scotch:BOOL=ON.");
  }
};

#endif // HAVE_ZOLTAN2_SCOTCH

////////////////////////////////////////////////////////////////////////

} // namespace Zoltan2

#endif
