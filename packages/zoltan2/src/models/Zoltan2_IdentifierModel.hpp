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

/*! \file Zoltan2_IdentifierModel.hpp
    \brief Defines the IdentifierModel interface.
*/

#ifndef _ZOLTAN2_IDENTIFIERMODEL_HPP_
#define _ZOLTAN2_IDENTIFIERMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_Adapter.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief IdentifierModel defines the interface for all identifier models.

    The constructor of the IdentifierModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.
*/

template <typename Adapter>
class IdentifierModel : public Model<Adapter> 
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef typename Adapter::gid_t     gid_t;
  typedef IdentifierMap<typename Adapter::user_t> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;
#endif

  /*! \brief Constructor
       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param modelFlags   bit map of Zoltan2::IdentifierModelFlags
   */
  
  IdentifierModel(const Adapter *ia, const RCP<const Environment> &env,
                  const RCP<const Comm<int> > &comm, modelFlag_t &modelFlags);

  /*! Returns the number of identifiers on this process.
   */
  inline size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  inline global_size_t getGlobalNumIdentifiers() const {
    return numGlobalIdentifiers_;
  }

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  inline int getIdentifierWeightDim() const { return this->getNumWeights(); }

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.
       \return The number of ids in the Ids list.
   */
  inline size_t getIdentifierList(ArrayView<const gno_t> &Ids,
                                  ArrayView<input_t> &wgts) const 
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, userWeightDim_);
    size_t n = getLocalNumIdentifiers();
    if (n){
      if (gnosAreGids_) Ids = gids_(0, n);
      else              Ids = gnosConst_(0, n);
    }
    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  inline size_t getLocalNumObjects() const {return getLocalNumIdentifiers();}

  inline size_t getGlobalNumObjects() const {return getGlobalNumIdentifiers();}

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  int userWeightDim_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

////////////////////////////////////////////////////
template <typename Adapter>
  IdentifierModel<Adapter>::IdentifierModel( 
    const Adapter *ia,
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), userWeightDim_(0), weights_(), gnos_(), gnosConst_()
{
  // Get the local and global problem size
  size_t nLocalIds = ia->getLocalNumIDs();
  gno_t lsum = nLocalIds;
  reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum,
    &numGlobalIdentifiers_);

  // Get the number of weights
  // Use max weight dim over all processes as userWeightDim_
  int tmp = ia->getNumWeightsPerID();
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 1,
      &tmp, &userWeightDim_);

  // Prepare to store views from input adapter
  // TODO:  Do we have to store these views, or can we get them on an 
  // TODO:  as-needed basis?
  Array<const scalar_t *> wgts(userWeightDim_, (const scalar_t *)NULL);
  Array<int> wgtStrides(userWeightDim_, 0);
  Array<lno_t> weightArrayLengths(userWeightDim_, 0);

  if (userWeightDim_ > 0){
    input_t *w = new input_t [userWeightDim_];
    weights_ = arcp<input_t>(w, 0, userWeightDim_);
  }

  const gid_t *gids=NULL;
  
  // Get the input adapter's views
  try{
    ia->getIDsView(gids);
    for (int idx=0; idx < userWeightDim_; idx++)
      ia->getWeightsView(wgts[idx], wgtStrides[idx], idx);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (nLocalIds){
    gids_ = arcp(gids, 0, nLocalIds, false);

    if (userWeightDim_ > 0){
      for (int i=0; i < userWeightDim_; i++){
        if (wgts[i] != NULL){
          ArrayRCP<const scalar_t> wgtArray(
            wgts[i], 0, nLocalIds*wgtStrides[i], false);
          weights_[i] = input_t(wgtArray, wgtStrides[i]);
          weightArrayLengths[i] = nLocalIds;
        }
      }
    }
  }

  this->setWeightArrayLengths(weightArrayLengths, *comm_);

  // TODO:  Why does an IdentifierModel need an IdentifierMap?
  // TODO:  Currently is useful only if gid_t is not Teuchos::Ordinal
  RCP<const idmap_t> idMap;
  try{
    if (modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE))
      idMap = rcp(new idmap_t(env_, comm_, gids_, true));
    else
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
  }
  Z2_FORWARD_EXCEPTIONS;

  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

  if (!gnosAreGids_ && nLocalIds>0){
    gno_t *tmpGno = new gno_t [nLocalIds];
    env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, tmpGno);
    gnos_ = arcp(tmpGno, 0, nLocalIds);

    try{
      ArrayRCP<gid_t> gidsNonConst = arcp_const_cast<gid_t>(gids_);
      idMap->gidTranslate( gidsNonConst(0,nLocalIds),  gnos_(0,nLocalIds),
        TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);

  env_->memory("After construction of identifier model");
}

}  // namespace Zoltan2

#endif
