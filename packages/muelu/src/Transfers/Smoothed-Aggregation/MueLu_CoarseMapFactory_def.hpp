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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_CoarseMapFactory_def.hpp
 *
 *  Created on: Oct 12, 2012
 *      Author: wiesner
 */

#ifndef MUELU_COARSEMAPFACTORY_DEF_HPP_
#define MUELU_COARSEMAPFACTORY_DEF_HPP_

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_CoarseMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoarseMapFactory(RCP<const FactoryBase> aggregatesFact, RCP<const FactoryBase> nullspaceFact)
  : aggregatesFact_(aggregatesFact), nullspaceFact_(nullspaceFact), domainGidOffset_(0)
  {
    stridedBlockId_ = -1; // default: blocked map with constant blocksize "NSDim"
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~CoarseMapFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
    currentLevel.DeclareInput("Nullspace",  nullspaceFact_.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "CoarseMap factory", currentLevel);

    //TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::NullspaceFactory::Build(): NullspaceFactory can be used for finest level (LevelID == 0) only.");

    RCP<Aggregates>  aggregates = currentLevel.Get< RCP<Aggregates> >("Aggregates", aggregatesFact_.get());
    GlobalOrdinal numAggs = aggregates->GetNumAggregates();

    // extract communicator
    RCP<const Teuchos::Comm<int> > comm = aggregates->GetMap()->getComm();

    // determine nullspace dimension
    RCP<MultiVector> nullspace  = currentLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());
    const size_t NSDim = nullspace->getNumVectors();

    // check for consistency of striding information with NSDim and nCoarseDofs
    if( stridedBlockId_== -1 ) {
      // this means we have no real strided map but only a block map with constant blockSize "NSDim"
      TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() > 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");
      stridingInfo_.clear();
      stridingInfo_.push_back(NSDim);
      TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() != 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");
    } else {
      // stridedBlockId_ > -1, set by user
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId_ > Teuchos::as<LO>(stridingInfo_.size() - 1) , Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): it is stridingInfo_.size() <= stridedBlockId_. error.");
      size_t stridedBlockSize = stridingInfo_[stridedBlockId_];
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockSize != NSDim , Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): dimension of strided block != NSDim. error.");
    }

    GetOStream(Statistics1, 0) << "domainGIDOffset: " << domainGidOffset_ << " block size: " << getFixedBlockSize() << " stridedBlockId: " << stridedBlockId_ << std::endl;

    // number of coarse level dofs (fixed by number of aggregates and blocksize data)
    GlobalOrdinal nCoarseDofs = numAggs * getFixedBlockSize();
    GlobalOrdinal indexBase   = aggregates->GetMap()->getIndexBase();

    RCP<const Map> coarseMap = StridedMapFactory::Build(aggregates->GetMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        nCoarseDofs,
        indexBase,
        stridingInfo_,
        comm,
        stridedBlockId_,
        domainGidOffset_);

    currentLevel.Set("CoarseMap", coarseMap, this);
  } // Build

} //namespace MueLu

#endif /* MUELU_COARSEMAPFACTORY_DEF_HPP_ */
