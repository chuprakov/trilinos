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
#ifndef MUELU_AMALGAMATIONFACTORY_DEF_HPP
#define MUELU_AMALGAMATIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_AmalgamationFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A"); // sub-block from blocked A
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
  {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    LO fullblocksize    = 1;    // block dim for fixed size blocks
    GO offset           = 0;   // global offset of dof gids
    LO blockid          = -1;  // block id in strided map
    LO nStridedOffset   = 0;   // DOF offset for strided block id "blockid" (default = 0)
    LO stridedblocksize = fullblocksize; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)
    GO indexBase        = A->getRowMap()->getIndexBase();  // index base for maps

    // 1) check for blocking/striding information

    if (A->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // NOTE: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      RCP<const StridedMap> stridedRowMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
      TEUCHOS_TEST_FOR_EXCEPTION(stridedRowMap == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
      fullblocksize = stridedRowMap->getFixedBlockSize();
      offset        = stridedRowMap->getOffset();
      blockid       = stridedRowMap->getStridedBlockId();

      if (blockid > -1) {
        std::vector<size_t> stridingInfo = stridedRowMap->getStridingData();
        for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
          nStridedOffset += stridingInfo[j];
        stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);

      } else {
        stridedblocksize = fullblocksize;
      }
      oldView = A->SwitchToView(oldView);
      GetOStream(Runtime1, 0) << "AmalagamationFactory::Build():" << " found fullblocksize=" << fullblocksize << " and stridedblocksize=" << stridedblocksize << " from strided maps. offset=" << offset << std::endl;

    } else {
      GetOStream(Warnings0, 0) << "AmalagamationFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;
    }
    // TODO: maybe no striding information on coarser levels -> misuse nullspace vector?

    // 2) prepare maps for amalgamated graph of A and
    //    setup unamalgamation information

    RCP<std::vector<GlobalOrdinal> > gNodeIds; // contains global node ids on current proc
    gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    gNodeIds->empty();

    // in nodegid2dofgids_ for each node on the current proc a vector of
    // the corresponding DOFs gids is stored.
    // The map contains all nodes the current proc has connections to (including
    // nodes that are stored on other procs when there are off-diagonal entries in A)
    nodegid2dofgids_ = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);

    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = A->getColMap();

    GlobalOrdinal cnt_amalRows = 0; // counts number of nodes (rows in amalgamated matrix) on current proc
    LocalOrdinal nColEle = Teuchos::as<LocalOrdinal>(colMap->getNodeNumElements());
    for (LocalOrdinal i = 0; i < nColEle; i++) {
      // get global DOF id
      GlobalOrdinal gDofId = colMap->getGlobalElement(i);

      // translate DOFGid to node id
      GlobalOrdinal gNodeId = DOFGid2NodeId(gDofId, A, fullblocksize, offset, indexBase);

      // gblockid -> gDofId/lDofId
      if (nodegid2dofgids_->count(gNodeId) == 0) {

        // current column DOF gDofId belongs to a node that has not been added
        // to nodeid2dofgids_ yet. Do it now and add ALL DOFs of node gNodeId to
        // unamalgamation information.
        // Note: we use offset and fullblocksize, ie. information from strided maps indirectly
        std::vector<GlobalOrdinal> DOFs;

        DOFs.reserve(stridedblocksize);
        for (LocalOrdinal k = 0; k < stridedblocksize; k++) {
          // here, the assumption is, that the node map has the same indexBase as the dof map
          //                            this is the node map index base                    this is the dof map index base
          GO gDofIndex = offset + (gNodeId-indexBase)*fullblocksize + nStridedOffset + k + indexBase;
          if (colMap->isNodeGlobalElement(gDofIndex))
            DOFs.push_back(gDofIndex);
        }

        (*nodegid2dofgids_)[gNodeId] = DOFs;

        if (rowMap->isNodeGlobalElement(gDofId)) {
          gNodeIds->push_back(gNodeId);
          cnt_amalRows++; // new local block row in amalgamated matrix graph
        }
      }
    }

    // store (un)amalgamation information on current level
    RCP<AmalgamationInfo> amalgamationData = rcp(new AmalgamationInfo());
    amalgamationData->SetAmalgamationParams(nodegid2dofgids_);
    amalgamationData->SetNodeGIDVector(gNodeIds);
    amalgamationData->SetNumberOfNodes(cnt_amalRows);
    Set(currentLevel, "UnAmalgamationInfo", amalgamationData);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const GlobalOrdinal AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DOFGid2NodeId(GlobalOrdinal gid, const RCP<Matrix>& A, LocalOrdinal blockSize, const GlobalOrdinal offset, const GlobalOrdinal indexBase) {
    // here, the assumption is, that the node map has the same indexBase as the dof map
    GlobalOrdinal globalblockid = ((GlobalOrdinal) gid - offset - indexBase) / blockSize + indexBase;
    return globalblockid;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UnamalgamateAggregates(const Aggregates& aggregates,
const AmalgamationInfo& amalgInfo, Teuchos::ArrayRCP<LocalOrdinal> & aggStart, Teuchos::ArrayRCP<GlobalOrdinal> & aggToRowMap) {
    int myPid = aggregates.GetMap()->getComm()->getRank();
    //const Map& map = *(aggregates.GetMap());
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();

    GO total=0;
    std::vector<LO> sizes(aggregates.GetNumAggregates());
    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];
      if (procWinner[lnode] == myPid) {
        //GO gnodeid = map.getGlobalElement(lnode);
        GO gnodeid = (aggregates.GetMap())->getGlobalElement(lnode);
        std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
        total += Teuchos::as<LO>(gDofIds.size());
        sizes[myAgg] += Teuchos::as<LO>(gDofIds.size());
      }
    }
    aggToRowMap = ArrayRCP<GO>(total,0);

    aggStart = ArrayRCP<LO>(aggregates.GetNumAggregates()+1,0);
    aggStart[0]=0;
    for (GO i=0; i<aggregates.GetNumAggregates(); ++i) {
      aggStart[i+1] = aggStart[i] + sizes[i];
    }

    // count, how many dofs have been recorded for each aggregate so far
    Array<LO> numDofs(aggregates.GetNumAggregates(), 0); // empty array with number of Dofs for each aggregate

    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];

      if (procWinner[lnode] == myPid) {
        //GO gnodeid = map.getGlobalElement(lnode);
        GO gnodeid = (aggregates.GetMap())->getGlobalElement(lnode);
        std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
        LO gDofIds_size = Teuchos::as<LO>(gDofIds.size());
        for (LO gDofId=0; gDofId < gDofIds_size; ++gDofId) {
          //aggToRowMap[ aggStart[myAgg] + gDofId ] = gDofIds[gDofId]; // fill aggToRowMap structure
          aggToRowMap[ aggStart[myAgg] + numDofs[myAgg] ] = gDofIds[gDofId]; // fill aggToRowMap structure
          ++(numDofs[myAgg]);
        }
      }
    }
    // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

  } //UnamalgamateAggregates

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeUnamalgamatedImportDofMap(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo) {
    Teuchos::RCP<const Map> nodeMap = aggregates.GetMap(); //aggregates.GetVertex2AggId();

    Teuchos::RCP<std::vector<GO> > myDofGids = Teuchos::rcp(new std::vector<GO>);
    LO nodeElements = Teuchos::as<LO>(nodeMap->getNodeNumElements());
    for(LO n = 0; n<nodeElements; n++) {
      GO gnodeid = (GO) nodeMap->getGlobalElement(n);
      std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
      for(typename std::vector<GO>::iterator gDofIdsIt = gDofIds.begin(); gDofIdsIt != gDofIds.end(); gDofIdsIt++) {
        myDofGids->push_back(*gDofIdsIt);
      }
    }

    Teuchos::ArrayRCP<GO> arr_myDofGids = Teuchos::arcp( myDofGids );
    Teuchos::RCP<Map> importDofMap = MapFactory::Build(aggregates.GetMap()->lib(),
                                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), arr_myDofGids(),
                                                       aggregates.GetMap()->getIndexBase(), aggregates.GetMap()->getComm());
    return importDofMap;
  }

} //namespace MueLu

#endif /* MUELU_SUBBLOCKUNAMALGAMATIONFACTORY_DEF_HPP */

