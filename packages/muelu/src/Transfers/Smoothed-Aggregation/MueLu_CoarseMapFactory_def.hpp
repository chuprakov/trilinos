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
 * MueLu_CoarseMapFactory_def.hpp
 *
 *  Created on: Oct 12, 2012
 *      Author: wiesner
 */

#ifndef MUELU_COARSEMAPFACTORY_DEF_HPP_
#define MUELU_COARSEMAPFACTORY_DEF_HPP_

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include <Teuchos_Array.hpp>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_CoarseMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoarseMapFactory()
  : domainGidOffset_(0)
  {
    //stridedBlockId_ = -1; // default: blocked map with constant blocksize "NSDim"
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~CoarseMapFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Generating factory for aggregates.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",  Teuchos::null, "Generating factory for null space.");

    validParamList->set< std::string  >("Striding info", "{}", "Striding information");
    validParamList->set< LocalOrdinal >("Strided block id", -1, "Strided block id");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "Aggregates");
    Input(currentLevel, "Nullspace");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::setDomainMapOffset(GlobalOrdinal offset) {
    TEUCHOS_TEST_FOR_EXCEPTION(offset < 0, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::SetDomainMapOffset: domain map offset for coarse gids of tentative prolongator must not be smaller than zero. Error.");
    domainGidOffset_ = offset;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getDomainMapOffset() const {
    return domainGidOffset_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::setStridingData(std::vector<size_t> stridingInfo) {
    // store striding map in internal variable
    stridingInfo_ = stridingInfo;

    // try to remove string "Striding info" from parameter list to make sure,
    // that stridingInfo_ is not replaced in the Build call.
    std::string strStridingInfo; strStridingInfo.clear();
    SetParameter("Striding info", ParameterEntry(strStridingInfo));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Aggregates>  aggregates = Get< RCP<Aggregates> >(currentLevel, "Aggregates");
    RCP<MultiVector> nullspace  = Get< RCP<MultiVector> >(currentLevel, "Nullspace");

    GlobalOrdinal                  numAggs = aggregates->GetNumAggregates();
    const size_t                   NSDim   = nullspace->getNumVectors();
    RCP<const Teuchos::Comm<int> > comm    = aggregates->GetMap()->getComm();

    // read in stridingInfo from parameter list and fill the internal member variable
    // read the data only if the parameter "Striding info" exists and is non-empty
    const ParameterList & pL = GetParameterList();
    if(pL.isParameter("Striding info")) {
      std::string strStridingInfo = pL.get<std::string>("Striding info");
      if(strStridingInfo.empty() == false) {
        Teuchos::Array<size_t> arrayVal = Teuchos::fromStringToArray<size_t>(strStridingInfo);
        stridingInfo_ = Teuchos::createVector(arrayVal);
      }
    }

    LocalOrdinal stridedBlockId = pL.get<LocalOrdinal>("Strided block id");

    // check for consistency of striding information with NSDim and nCoarseDofs
    if (stridedBlockId== -1) {
      // this means we have no real strided map but only a block map with constant blockSize "NSDim"
      TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() > 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");
      stridingInfo_.clear();
      stridingInfo_.push_back(NSDim);
      TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() != 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");

    } else {
      // stridedBlockId > -1, set by user
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId > Teuchos::as<LO>(stridingInfo_.size() - 1) , Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): it is stridingInfo_.size() <= stridedBlockId_. error.");
      size_t stridedBlockSize = stridingInfo_[stridedBlockId];
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockSize != NSDim , Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): dimension of strided block != NSDim. error.");
    }

    GetOStream(Statistics2) << "domainGIDOffset: " << domainGidOffset_ << " block size: " << getFixedBlockSize() << " stridedBlockId: " << stridedBlockId << std::endl;

    // number of coarse level dofs (fixed by number of aggregates and blocksize data)
    GlobalOrdinal nCoarseDofs = numAggs * getFixedBlockSize();
    GlobalOrdinal indexBase   = aggregates->GetMap()->getIndexBase();

    RCP<const Map> coarseMap = StridedMapFactory::Build(aggregates->GetMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        nCoarseDofs,
        indexBase,
        stridingInfo_,
        comm,
        stridedBlockId,
        domainGidOffset_);

    Set(currentLevel, "CoarseMap", coarseMap);
  } // Build

} //namespace MueLu

#endif /* MUELU_COARSEMAPFACTORY_DEF_HPP_ */
