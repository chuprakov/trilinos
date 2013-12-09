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

/*! \file Zoltan2_XpetraMultiVectorAdapter.hpp
    \brief Defines the XpetraMultiVectorAdapter 
*/

#ifndef _ZOLTAN2_XPETRAMULTIVECTORADAPTER_HPP_
#define _ZOLTAN2_XPETRAMULTIVECTORADAPTER_HPP_

#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_Util.hpp>

#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

namespace Zoltan2 {

/*!  \brief An adapter for Xpetra::MultiVector.

    The template parameter is the user's input object: 
    \li \c Epetra_MultiVector
    \li \c Tpetra::MultiVector
    \li \c Xpetra::MultiVector

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.
*/

template <typename User>
  class XpetraMultiVectorAdapter : public VectorAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef VectorAdapter<User>       base_adapter_t;
  typedef User user_t;

  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> x_mvector_t;
  typedef Xpetra::TpetraMultiVector<
    scalar_t, lno_t, gno_t, node_t> xt_mvector_t;
  typedef Xpetra::EpetraMultiVector xe_mvector_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraMultiVectorAdapter() { }

  /*! \brief Constructor   
   *
   *  \param invector  the user's Xpetra, Tpetra or Epetra MultiVector object
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per multivector element is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight dimension \c n for multivector element
   *     \c k should be found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  XpetraMultiVectorAdapter(const RCP<const User> &invector,
    vector<const scalar_t *> &weights, vector<int> &weightStrides);

  /*! \brief Constructor for case when weights are not being used.
   *
   *  \param invector  the user's Xpetra, Tpetra or Epetra MultiVector object
   */

  XpetraMultiVectorAdapter(const RCP<const User> &invector);


  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  size_t getLocalNum() const { return vector_->getLocalLength();}

  size_t getIDsView(const gid_t *&ids) const
  { 
    ids = map_->getNodeElementList().getRawPtr();
    return vector_->getLocalLength();
  }

  int getNumWeightsPer() const { return numWeights_;}

  size_t getWeightsView(const scalar_t *&weights, int &stride, int idx) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight index",
      idx >= 0 && idx < numWeights_, BASIC_ASSERTION);
    size_t length;
    weights_[idx].getStridedList(length, weights, stride);
    return length;
  }

  ////////////////////////////////////////////////////
  // The VectorAdapter interface.
  ////////////////////////////////////////////////////

  int getNumVectors() const {return vector_->getNumVectors();}

  size_t getVectorView(const scalar_t *&elements, int &stride, int idx=0) const;

  template <typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<const User> invector_;
  RCP<const x_mvector_t> vector_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map_;
  RCP<Environment> env_;    // for error messages, etc.
  lno_t base_;

  int numWeights_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;
};

//////////////////////////////////////////////////////////
// Definitions
//////////////////////////////////////////////////////////

template <typename User>
  XpetraMultiVectorAdapter<User>::XpetraMultiVectorAdapter(
    const RCP<const User> &invector,
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      invector_(invector), vector_(), map_(), 
      env_(rcp(new Environment)), base_(),
      numWeights_(weights.size()), weights_(weights.size())
{
  typedef StridedData<lno_t, scalar_t> input_t;

  vector_ = XpetraTraits<User>::convertToXpetra(invector);
  map_ = vector_->getMap();
  base_ = map_->getIndexBase();

  size_t length = vector_->getLocalLength();

  if (length > 0 && numWeights_ > 0){
    int stride = 1;
    for (int w=0; w < numWeights_; w++){
      if (weightStrides.size())
        stride = weightStrides[w];
      ArrayRCP<const scalar_t> wgtV(weights[w], 0, stride*length, false); 
      weights_[w] = input_t(wgtV, stride);
    }
  }
}


template <typename User>
  XpetraMultiVectorAdapter<User>::XpetraMultiVectorAdapter(
    const RCP<const User> &invector):
      invector_(invector), vector_(), map_(), 
      env_(rcp(new Environment)), base_(),
      numWeights_(0), weights_()
{
  vector_ = XpetraTraits<User>::convertToXpetra(invector);
  map_ = vector_->getMap();
  base_ = map_->getIndexBase();
}

template <typename User>
  size_t XpetraMultiVectorAdapter<User>::getVectorView(
    const scalar_t *&elements, int &stride, int idx) const
{
  size_t vecsize;
  stride = 1;
  elements = NULL;
  if (map_->lib() == Xpetra::UseTpetra){
    const xt_mvector_t *tvector = 
      dynamic_cast<const xt_mvector_t *>(vector_.get());
     
    vecsize = tvector->getLocalLength();
    if (vecsize > 0){
      ArrayRCP<const scalar_t> data = tvector->getData(idx);
      elements = data.get();
    }
  }
  else if (map_->lib() == Xpetra::UseEpetra){
    const xe_mvector_t *evector = 
      dynamic_cast<const xe_mvector_t *>(vector_.get());
      
    vecsize = evector->getLocalLength();
    if (vecsize > 0){
      ArrayRCP<const double> data = evector->getData(idx);

      // Cast so this will compile when scalar_t is not double,
      // a case when this code should never execute.
      elements = reinterpret_cast<const scalar_t *>(data.get());
    }
  }
  else{
    throw logic_error("invalid underlying lib");
  }

  ArrayView<const gid_t> gids = map_->getNodeElementList();
  return vecsize;
}

template <typename User>
  template <typename Adapter>
    size_t XpetraMultiVectorAdapter<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<Adapter> &solution) const
{
  size_t len = solution.getLocalNumberOfIds();
  const gid_t *gids = solution.getIdList();
  const partId_t *parts = solution.getPartList();
  ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false);
  ArrayRCP<partId_t> partList = arcp(const_cast<partId_t *>(parts), 0, len, 
    false);
  ArrayRCP<lno_t> dummyIn;
  ArrayRCP<gid_t> importList;
  ArrayRCP<lno_t> dummyOut;
  size_t numNewRows;

  const RCP<const Comm<int> > comm = map_->getComm();

  try{
    // Get an import list
    numNewRows = solution.convertSolutionToImportList(
      0, dummyIn, importList, dummyOut);
  }
  Z2_FORWARD_EXCEPTIONS;

  RCP<const User> inPtr = rcp(&in, false);
  lno_t localNumElts = numNewRows;

  RCP<const User> outPtr = XpetraTraits<User>::doMigration(
   inPtr, localNumElts, importList.get());

  out = const_cast<User *>(outPtr.get());
  outPtr.release();
  return numNewRows;
}
  
}  //namespace Zoltan2
  
#endif
