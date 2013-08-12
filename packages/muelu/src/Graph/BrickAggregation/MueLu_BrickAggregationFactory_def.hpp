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
#ifndef MUELU_BRICKAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_BRICKAGGREGATIONFACTORY_DEF_HPP_

#include "MueLu_BrickAggregationFactory_decl.hpp"
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>


#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for coordinates");
    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory for matrix");
    validParamList->set< int >                   ("bx",                             2, "Number of brick points for x axis");
    validParamList->set< int >                   ("by",                             2, "Number of brick points for x axis");
    validParamList->set< int >                   ("bz",                             2, "Number of brick points for x axis");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");
  }

  // The current implementation cannot deal with bricks larger than 3x3(x3).
  // The reason is that aggregation infrastructure in place has major drawbacks.
  //
  // Aggregates class is constructed with a help of a provided map, either taken from
  // a graph, or provided directly. This map is usually taken to be a column map of a
  // matrix. The reason for that is that if we have an overlapped aggregation, we want
  // the processor owning aggregates to store agg id for all nodes in this aggregate.
  // If we used row map, there would be no way for the processor to know whether there
  // are some other nodes on a different processor which belong to its aggregate.
  // On the other hand, using column map allows both vertex2AggId and procWinner arrays
  // in Aggregates class to store some extra data, such as whether nodes belonging to a
  // different processor belong to this processor aggregate.
  // The drawback of this is that it stores only overlap=1 data. For aggressive coarsening,
  // such a brick aggregation with a large single dimension of brick, it could happen
  // that we need to know depth two or more extra nodes in the other processor subdomain.
  //
  // Another issue is that we may have some implicit connection between aggregate map and
  // maps of A used in the construction of a tentative prolongator.
  //
  // Another issue is that it seems that some info is unused or not required. Specifically,
  // it seems that if a node belongs to an aggregate on a different processor, we don't
  // actually need to set vertex2AggId and procWinner, despite the following comment in
  // Aggregates decl:
  //      vertex2AggId[k] gives a local id
  //      corresponding to the aggregate to which
  //      local id k has been assigned.  While k
  //      is the local id on my processor (MyPID)
  //      vertex2AggId[k] is the local id on the
  //      processor which actually owns the
  //      aggregate. This owning processor has id
  //      given by procWinner[k].
  // It is possible that that info is only used during arbitration in CoupledAggregationFactory.
  //
  // The steps that we need to do to resolve this issue:
  //   - Break the link between maps in TentativePFactory, allowing any maps in Aggregates
  //   - Allow Aggregates to construct their own maps, if necessary, OR
  //   - construct aggregates based on row map
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    const ParameterList& pL = GetParameterList();

    RCP<MultiVector> coords = Get< RCP<MultiVector> >(currentLevel, "Coordinates");
    RCP<Matrix>      A      = Get< RCP<Matrix> >     (currentLevel, "A");
    RCP<const Map>   rowMap = A->getRowMap();
    RCP<const Map>   colMap = A->getColMap();

    RCP<const Teuchos::Comm<int> > comm = rowMap->getComm();
    int myRank = comm->getRank();

    int numPoints = colMap->getNodeNumElements();

    bx_ = pL.get<int>("bx");
    by_ = pL.get<int>("by");
    bz_ = pL.get<int>("bz");

    // TODO: deal with block size > 1
    TEUCHOS_TEST_FOR_EXCEPTION(bx_ > 3 || by_ > 3 || bz_ > 3, Exceptions::RuntimeError, "Currently cannot deal with brick size > 3");

    RCP<MultiVector> overlappedCoords = coords;
    RCP<const Import> importer = ImportFactory::Build(coords->getMap(), colMap);
    if (!importer.is_null()) {
      overlappedCoords = MultiVectorFactory::Build(colMap, coords->getNumVectors());
      overlappedCoords->doImport(*coords, *importer, Xpetra::INSERT);
    }

    // Setup misc structures
    // Logically, we construct enough data to query topological information of a rectangular grid
    Setup(comm, overlappedCoords, colMap);

    GetOStream(Runtime0,0) << "Using brick size: " << bx_ << " x " << by_ << (nDim_ == 3 ? "x " + toString(bz_) : "") << std::endl;

    // Construct aggregates
    RCP<Aggregates> aggregates = rcp(new Aggregates(colMap));
    aggregates->setObjectLabel("Brick");

    Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
    Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates->GetProcWinner()  ->getDataNonConst(0);

    // In the first pass, we set a mapping from a vertex to aggregate global id. We deal with a structured
    // rectangular mesh, therefore we know the structure of aggregates. For each vertex we can tell exactly
    // which aggregate it belongs to.
    // If we determine that the aggregate does not belong to us (i.e. the root vertex does not belong to this
    // processor, or is outside and we lost "" arbitration), we record the global aggregate id in order to
    // fetch the local info from the processor owning the aggregate. This is required for aggregates, as it
    // uses the local aggregate ids of the owning processor.
    std::set<GlobalOrdinal> myAggGIDs, remoteAggGIDs;
    for (LocalOrdinal LID = 0; LID < numPoints; LID++) {
      GlobalOrdinal aggGID = getAggGID(LID);

      if ((revMap_.find(getRoot(LID)) != revMap_.end()) && rowMap->isNodeGlobalElement(colMap->getGlobalElement(revMap_[getRoot(LID)]))) {
        // Root of the brick aggregate containing GID (<- LID) belongs to us
        vertex2AggId[LID] = aggGID;
        myAggGIDs.insert(aggGID);

        if (isRoot(LID))
          aggregates->SetIsRoot(LID);

      } else {
        remoteAggGIDs.insert(aggGID);
      }
    }
    size_t numAggregates = myAggGIDs    .size();
    size_t numRemote     = remoteAggGIDs.size();
    aggregates->SetNumAggregates(numAggregates);

    std::map<GlobalOrdinal,LocalOrdinal> AggG2L;  // Map: Agg GID -> Agg LID (possibly on a different processor)
    std::map<GlobalOrdinal,int>          AggG2R;  // Map: Agg GID -> processor rank owning aggregate

    Array<GlobalOrdinal> myAggGIDsArray(numAggregates), remoteAggGIDsArray(numRemote);

    // Fill in the maps for aggregates that we own
    size_t ind = 0;
    for (typename std::set<GlobalOrdinal>::const_iterator it = myAggGIDs.begin(); it != myAggGIDs.end(); it++) {
      AggG2L[*it] = ind;
      AggG2R[*it] = myRank;

      myAggGIDsArray[ind++] = *it;
    }

    // The map is a convenient way to fetch remote local indices from global indices.
    RCP<Map> aggMap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                                 myAggGIDsArray, 0, comm);

    ind = 0;
    for (typename std::set<GlobalOrdinal>::const_iterator it = remoteAggGIDs.begin(); it != remoteAggGIDs.end(); it++)
      remoteAggGIDsArray[ind++] = *it;

    // Fetch the required aggregate local ids and ranks
    Array<int>          remoteProcIDs(numRemote);
    Array<LocalOrdinal> remoteLIDs   (numRemote);
    aggMap->getRemoteIndexList(remoteAggGIDsArray, remoteProcIDs, remoteLIDs);

    // Fill in the maps for aggregates that we don't own but which have some of our vertices
    for (size_t i = 0; i < numRemote; i++) {
      AggG2L[remoteAggGIDsArray[i]] = remoteLIDs   [i];
      AggG2R[remoteAggGIDsArray[i]] = remoteProcIDs[i];
    }

    // Remap aggregate GIDs to LIDs and set up owning processors
    for (LocalOrdinal LID = 0; LID < numPoints; LID++) {
      if (revMap_.find(getRoot(LID)) != revMap_.end() && rowMap->isNodeGlobalElement(colMap->getGlobalElement(revMap_[getRoot(LID)]))) {
        GlobalOrdinal aggGID = vertex2AggId[LID];

        vertex2AggId[LID] = AggG2L[aggGID];
        procWinner  [LID] = AggG2R[aggGID];
      }
    }

    GlobalOrdinal numGlobalRemote;
    sumAll(comm, Teuchos::as<GlobalOrdinal>(numRemote), numGlobalRemote);
    aggregates->AggregatesCrossProcessors(numGlobalRemote);

    Set(currentLevel, "Aggregates", aggregates);

    GetOStream(Statistics0, 0) << aggregates->description() << std::endl;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(const RCP<const Teuchos::Comm<int> >& comm, const RCP<MultiVector>& coords, const RCP<const Map>& map) const {
    nDim_ = coords->getNumVectors();

    x_    = coords->getData(0);
    xMap_ = Construct1DMap(comm, x_);
    nx_   = xMap_->size();

    y_    = coords->getData(1);
    yMap_ = Construct1DMap(comm, y_);
    ny_   = yMap_->size();

    nz_   = 1;
    if (nDim_ == 3) {
      z_    = coords->getData(2);
      zMap_ = Construct1DMap(comm, z_);
      nz_   = zMap_->size();
    }

    for (size_t ind = 0; ind < coords->getLocalLength(); ind++) {
      GlobalOrdinal i = (*xMap_)[(coords->getData(0))[ind]];
      GlobalOrdinal j = (*yMap_)[(coords->getData(1))[ind]];
      GlobalOrdinal k = 0;
      if (nDim_ == 3)
        k = (*zMap_)[(coords->getData(2))[ind]];

      revMap_[k*ny_*nx_ + j*nx_ + i] = ind;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<typename BrickAggregationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::container> BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Construct1DMap(const RCP<const Teuchos::Comm<int> >& comm, const ArrayRCP<const double>& x) const {
    int numProcs = comm->getSize();
    int n = x.size();

    RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    bool isMPI = !mpiComm.is_null();

    // Step 1: Create a local vector with unique coordinate points
    RCP<container> gMap = rcp(new container);
    for (int i = 0; i < n; i++)
      (*gMap)[x[i]] = 0;

    // Step 2: exchange coordinates
    Array<double>  recvBuf;
    int recvSize = 0;
    if (isMPI && numProcs > 1) {
      MPI_Comm rawComm;
      if (isMPI)
        rawComm = (*mpiComm->getRawMpiComm())();

      int           sendCnt = gMap->size(), cnt = 0;
      Array<int>    recvCnt(numProcs), Displs(numProcs);
      Array<double> sendBuf, recvBuf;

      sendBuf.resize(sendCnt);
      for (typename container::const_iterator cit = gMap->begin(); cit != gMap->end(); cit++)
        sendBuf[cnt++] = cit->first;

      MPI_Allgather(&sendCnt, 1, MPI_INT, recvCnt.getRawPtr(), 1, MPI_INT, rawComm);
      Displs[0] = 0;
      for (int i = 0; i < numProcs-1; i++)
        Displs[i+1] = Displs[i] + recvCnt[i];
      recvSize = Displs[numProcs-1] + recvCnt[numProcs-1];
      recvBuf.resize(recvSize);
      MPI_Allgatherv(sendBuf.getRawPtr(), sendCnt, MPI_DOUBLE, recvBuf.getRawPtr(), recvCnt.getRawPtr(), Displs.getRawPtr(), MPI_DOUBLE, rawComm);

      for (int i = 0; i < recvSize; i++)
        (*gMap)[recvBuf[i]] = 0;
    }

    GlobalOrdinal cnt = 0;
    for (typename container::iterator it = gMap->begin(); it != gMap->end(); it++)
      it->second = cnt++;

    return gMap;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::isRoot(LocalOrdinal LID) const {
    int i = (*xMap_)[x_[LID]];
    int j = (*yMap_)[y_[LID]];
    int k = 0;
    if (nDim_ == 3)
      k = (*zMap_)[z_[LID]];

    return (k*ny_*nx_ + j*nx_ + i) == getRoot(LID);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getRoot(LocalOrdinal LID) const {
    int i = ((*xMap_)[x_[LID]]/bx_)*bx_ + (bx_-1)/2;
    int j = ((*yMap_)[y_[LID]]/by_)*by_ + (by_-1)/2;
    int k = 0;
    if (nDim_ == 3)
      k = ((*zMap_)[z_[LID]]/bz_)*bz_ + (bz_-1)/2;

    // Check if the actual root is outside of the domain. If it is, project the root to the domain
    if (i > nx_) i = nx_;
    if (j > ny_) j = ny_;
    if (k > nz_) k = nz_;

    return k*ny_*nx_ + j*nx_ + i;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getAggGID(LocalOrdinal LID) const {
    int naggx = nx_/bx_ + (nx_ % bx_ ? 1 : 0), i = (*xMap_)[x_[LID]]/bx_;
    int naggy = ny_/by_ + (ny_ % by_ ? 1 : 0), j = (*yMap_)[y_[LID]]/by_;
    int                                        k = 0;
    if (nDim_ == 3)
      k = (*zMap_)[z_[LID]]/bz_;

    return k*naggy*naggx + j*naggx + i;
  }


} //namespace MueLu

#endif /* MUELU_BRICKAGGREGATIONFACTORY_DEF_HPP_ */
