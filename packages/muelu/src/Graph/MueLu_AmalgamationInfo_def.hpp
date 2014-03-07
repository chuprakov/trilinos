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
 * MueLu_AmalgamationInfo_def.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_DEF_HPP_
#define MUELU_AMALGAMATIONINFO_DEF_HPP_

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "MueLu_AmalgamationInfo_decl.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UnamalgamateAggregates(const Aggregates& aggregates,
        Teuchos::ArrayRCP<LocalOrdinal>& aggStart, Teuchos::ArrayRCP<GlobalOrdinal>& aggToRowMap) const {
    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayView<const GO> nodeGlobalElts = aggregates.GetMap()->getNodeElementList();
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();
    GO numAggregates = aggregates.GetNumAggregates();

    std::vector<LO> sizes(numAggregates);
    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];
      if (procWinner[lnode] == myPid) {
        GO gnodeid = nodeGlobalElts[lnode];
        std::vector<GO> gDofIds = ComputeGlobalDOFs(gnodeid);
        sizes[myAgg] += Teuchos::as<LO>(gDofIds.size());
      }
    }
    aggStart = ArrayRCP<LO>(numAggregates+1,0);
    aggStart[0]=0;
    for (GO i=0; i<numAggregates; ++i) {
      aggStart[i+1] = aggStart[i] + sizes[i];
    }
    aggToRowMap = ArrayRCP<GO>(aggStart[numAggregates],0);

    // count, how many dofs have been recorded for each aggregate so far
    Array<LO> numDofs(numAggregates, 0); // empty array with number of Dofs for each aggregate

    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];

      if (procWinner[lnode] == myPid) {
        GO gnodeid = nodeGlobalElts[lnode];
        std::vector<GO> gDofIds = ComputeGlobalDOFs(gnodeid);
        LO gDofIds_size = Teuchos::as<LO>(gDofIds.size());
        for (LO gDofId=0; gDofId < gDofIds_size; ++gDofId) {
          aggToRowMap[ aggStart[myAgg] + numDofs[myAgg] ] = gDofIds[gDofId]; // fill aggToRowMap structure
          ++(numDofs[myAgg]);
        }
      }
    }
    // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

  } //UnamalgamateAggregates

  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeUnamalgamatedImportDofMap(const Aggregates& aggregates) const {
    Teuchos::RCP<const Map> nodeMap = aggregates.GetMap();

    Teuchos::RCP<std::vector<GO> > myDofGids = Teuchos::rcp(new std::vector<GO>);
    LO nodeElements = Teuchos::as<LO>(nodeMap->getNodeNumElements());
    for(LO n = 0; n<nodeElements; n++) {
      GO gnodeid = (GO) nodeMap->getGlobalElement(n);
      std::vector<GO> gDofIds = ComputeGlobalDOFs(gnodeid);
      for(typename std::vector<GO>::iterator gDofIdsIt = gDofIds.begin(); gDofIdsIt != gDofIds.end(); gDofIdsIt++) {
        myDofGids->push_back(*gDofIdsIt);
      }
    }

    Teuchos::ArrayRCP<GO> arr_myDofGids = Teuchos::arcp( myDofGids );
    Teuchos::RCP<Map> importDofMap = MapFactory::Build(aggregates.GetMap()->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), arr_myDofGids(), aggregates.GetMap()->getIndexBase(), aggregates.GetMap()->getComm());
    return importDofMap;
  }

  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::vector<GlobalOrdinal> AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeGlobalDOFs(GlobalOrdinal const gNodeID) const {
    //optimization:  preallocate DOFs and reuse it
    std::vector<GlobalOrdinal> DOFs;

    DOFs.reserve(stridedblocksize_);
    for (LocalOrdinal k = 0; k < stridedblocksize_; k++) {
    // here, the assumption is, that the node map has the same indexBase as the dof map
    //                            this is the node map index base                    this is the dof map index base
      GlobalOrdinal gDofIndex = offset_ + (gNodeID-indexBase_)*fullblocksize_ + nStridedOffset_ + k + indexBase_;
      if (columnMap_->isNodeGlobalElement(gDofIndex))
        DOFs.push_back(gDofIndex);
    }
    return DOFs;
  }

} //namespace


#endif /* MUELU_AMALGAMATIONINFO_DEF_HPP_ */
