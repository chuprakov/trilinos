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
#ifndef MUELU_REPARTITIONFACTORY_DEF_HPP
#define MUELU_REPARTITIONFACTORY_DEF_HPP

#include <algorithm>
#include <iostream>
#include <sstream>

#include "MueLu_RepartitionFactory_decl.hpp" // TMP JG NOTE: before other includes, otherwise I cannot test the fwd declaration in _def

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Hashtable.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <MueLu_CoupledAggregationCommHelper.hpp>

#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<int>        ("startLevel",                  1, "First level at which repartitioning can possibly occur. Repartitioning at finer levels is suppressed");
    validParamList->set<LO>         ("minRowsPerProcessor",      1000, "Minimum number of rows over all processes. If any process falls below this, repartitioning is initiated");
    validParamList->set<double>     ("nonzeroImbalance",          1.2, "Imbalance threshold, below which repartitioning is initiated. Imbalance is measured by "
                                                                       "ratio of maximum nonzeros over all processes to minimum number of nonzeros over all processes");

    validParamList->set<bool>       ("remapPartitions",         false, "Perform partition remapping to minimize data movement");
    validParamList->set<bool>       ("alwaysKeepProc0",          true, "Always keep processor 0 in subcommunicator");

    validParamList->set< RCP<const FactoryBase> >("A",         Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Partition", Teuchos::null, "Factory of the partition");

    // By default, the generating factory is 'this, not Teuchos::null, which means that neither user defined data nor FactoryManager
    // are used, which is unusual. That is because this class actually computes the number of partitions internally.
    //
    // However, one can specify the number of partitions to override this behaviour. In order to manually set the "number of partitions"
    // entry in the level by doing one of these alternatives:
    //    *) level.Set("numbers of partitions");
    //       myRepartitionFact.set< RCP<FactoryBase> >("number of partitions", Teuchos::null);
    //    *) level.Set("numbers of partitions", myRepartitionFact)
    // and doing appropriate requests
    RCP<const FactoryBase> rcpThis = rcpFromRef(*this);
    validParamList->set< RCP<const FactoryBase> >("number of partitions", rcpThis, "(advanced) Factory computing the number of partition. By default, an appropriate number of partition is computed internally.");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Partition");
    Input(currentLevel, "number of partitions");
  }

  template<class T> class MpiTypeTraits            { public: static MPI_Datatype getType(); };
  template<>        class MpiTypeTraits<long>      { public: static MPI_Datatype getType() { return MPI_LONG;      } };
  template<>        class MpiTypeTraits<int>       { public: static MPI_Datatype getType() { return MPI_INT;       } };
  template<>        class MpiTypeTraits<short>     { public: static MPI_Datatype getType() { return MPI_SHORT;     } };
#ifdef TEUCHOS_HAVE_LONG_LONG_INT
  template<>        class MpiTypeTraits<long long> { public: static MPI_Datatype getType() { return MPI_LONG_LONG; } };
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    const Teuchos::ParameterList & pL = GetParameterList();
    // Access parameters here to make sure that we set the parameter entry flag to "used" even in case of short-circuit evaluation.
    // TODO (JG): I don't really know if we want to do this.
    const int    startLevel          = pL.get<int>   ("startLevel");
    const LO     minRowsPerProcessor = pL.get<LO>    ("minRowsPerProcessor");
    const double nonzeroImbalance    = pL.get<double>("nonzeroImbalance");
    const bool   remapPartitions     = pL.get<bool>  ("remapPartitions");
    const bool   keepProc0           = pL.get<bool>  ("alwaysKeepProc0");

    // TODO: We only need a CrsGraph. This class does not have to be templated on Scalar types.
    RCP<Matrix>            A         = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<const Map>         rowMap    = A->getRowMap();
    GO                     indexBase = rowMap->getIndexBase();
    Xpetra::UnderlyingLib  lib       = rowMap->lib();

    RCP<const Teuchos::Comm<int> > origComm = rowMap->getComm();
    RCP<const Teuchos::Comm<int> > comm = origComm->duplicate();
    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

    // ======================================================================================================
    // Determine whether partitioning is needed
    // ======================================================================================================
    // We repartition if the following expression is true: !test1 && !test2 && (test3 || test4)
    //
    // NOTE: most tests include some global communication, which is why we currently only do tests until we make
    // a decision on whether to repartition. However, there is value in knowing how "close" we are to having to
    // rebalance an operator. So, it would probably be beneficial to do and report *all* tests.
    bool test1 = false, test2 = false, test3 = false, test4 = false;
    std::string msg1, msg2, msg3, msg4;

    // Test1: skip repartitioning if current level is less than the specified minimum level for repartitioning
    if (currentLevel.GetLevelID() < startLevel) {
      test1 = true;

      msg1 = "\n  current level = " + toString(currentLevel.GetLevelID()) + ", first level where repartitioning can happen is " + toString(startLevel);
    }

    // Test 2: check whether A is actually distributed, i.e. more than one processor owns part of A
    // TODO: this global communication can be avoided if we store the information with the matrix (it is known when matrix is created)
    // TODO: further improvements could be achieved when we use subcommunicator for the active set. Then we only need to check its size
    if (!test1) {
      int numActiveProcesses = 0;
      sumAll(comm, Teuchos::as<int>((A->getNodeNumRows() > 0) ? 1 : 0), numActiveProcesses);

      if (numActiveProcesses == 1)
        test2 = true;

      msg2 = "\n  # processes with rows = " + toString(numActiveProcesses);
    }

    // Test3: check whether number of rows on any processor satisfies the minimum number of rows requirement
    // NOTE: Test2 ensures that repartitionning is not done when there is only one processor (it may or may not satisfy Test3)
    if (!test1 && !test2 && minRowsPerProcessor > 0) {
      LO numMyRows = Teuchos::as<LO>(A->getNodeNumRows()), minNumRows, LOMAX = Teuchos::OrdinalTraits<LO>::max();
      LO haveFewRows = (numMyRows < minRowsPerProcessor ? 1 : 0), numWithFewRows = 0;
      sumAll(comm, haveFewRows, numWithFewRows);
      minAll(comm, (numMyRows > 0 ? numMyRows : LOMAX), minNumRows);

      // TODO: we could change it to repartition only if the number of processors with numRows < minNumRows is larger than some
      // percentage of the total number. This way, we won't repartition if 2 out of 1000 processors don't have enough elements.
      // I'm thinking maybe 20% threshold. To implement, simply add " && numWithFewRows < .2*numProcs" to the if statement.
      if (numWithFewRows > 0)
        test3 = true;

      msg3 = "\n  min # rows per proc = " + toString(minNumRows) + ", min allowable = " + toString(minRowsPerProcessor);
    }

    // Test4: check whether the balance in the number of nonzeros per processor is greater than threshold
    if (!test1 && !test2 && !test3) {
      GO minNnz, maxNnz, numMyNnz = Teuchos::as<GO>(A->getNodeNumEntries());
      maxAll(comm, numMyNnz,                           maxNnz);
      minAll(comm, (numMyNnz > 0 ? numMyNnz : maxNnz), minNnz); // min nnz over all active processors
      double imbalance = Teuchos::as<double>(maxNnz)/minNnz;

      if (imbalance > nonzeroImbalance)
        test4 = true;

      msg4 = "\n  nonzero imbalance = " + toString(imbalance) + ", max allowable = " + toString(nonzeroImbalance);
    }

    if (test1 || test2 || (!test3 && !test4)) {
      GetOStream(Statistics0, 0) << "Repartitioning?  NO:";
      if      (test1) GetOStream(Statistics0,0) << msg1        << std::endl;
      else if (test2) GetOStream(Statistics0,0) << msg2        << std::endl;
      else            GetOStream(Statistics0,0) << msg3 + msg4 << std::endl;

      Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
      return;
    }

    GetOStream(Statistics0,0) << "Repartitioning? YES:" << msg3 + msg4 << std::endl;

    // ======================================================================================================
    // Calculate number of partitions
    // ======================================================================================================
    // FIXME Quick way to figure out how many partitions there should be (same algorithm as ML)
    // FIXME Should take into account nnz? Perhaps only when user is using min #nnz per row threshold.
    GO numPartitions;
    if (IsAvailable(currentLevel, "number of partitions")) {
      numPartitions = Get<GO>(currentLevel, "number of partitions");
      GetOStream(Warnings0, 0) << "Using user-provided \"number of partitions\", the performance is unknown" << std::endl;

    } else {
      if (Teuchos::as<GO>(A->getGlobalNumRows()) < minRowsPerProcessor) {
        // System is too small, migrate it to a single processor
        numPartitions = 1;

      } else {
        // Make sure that each processor has approximately minRowsPerProcessor
        numPartitions = A->getGlobalNumRows() / minRowsPerProcessor;
      }
      numPartitions = std::min(numPartitions, Teuchos::as<GO>(numProcs));

      Set(currentLevel, "number of partitions", numPartitions);
    }
    GetOStream(Statistics0, 0) << "Number of partitions to use = " << numPartitions << std::endl;

    // ======================================================================================================
    // Construct decomposition vector
    // ======================================================================================================
    RCP<GOVector> decomposition;
    if (numPartitions == 1) {
      // Trivial case: decomposition is the trivial one, all zeros. We skip the call to Zoltan_Interface
      // (this is mostly done to avoid extra output messages, as even if we didn't skip there is a shortcut
      // in Zoltan[12]Interface).
      // TODO: We can probably skip more work in this case (like building all extra data structures)
      GetOStream(Warnings0, 0) << "Only one partition: Skip call to the repartitioner." << std::endl;
      decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), true);

    } else {
      decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");

      if (decomposition.is_null()) {
        GetOStream(Warnings0, 0) << "No repartitioning necessary: partitions were left unchanged by the repartitioner" << std::endl;
        Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
        return;
      }
    }

    // ======================================================================================================
    // Remap if necessary
    // ======================================================================================================
    if (remapPartitions) {
      SubFactoryMonitor m1(*this, "DeterminePartitionPlacement", currentLevel);

      DeterminePartitionPlacement(*A, *decomposition, numPartitions, keepProc0);
    }

    // ======================================================================================================
    // Construct importer
    // ======================================================================================================
    // At this point, the following is true:
    //  * Each processors owns 0 or 1 partitions
    //  * If a processor owns a partition, that partition number is equal to the processor rank
    //  * The decomposition vector contains the partitions ids that the corresponding GID belongs to

    ArrayRCP<const GO> decompEntries;
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);

#ifdef HAVE_MUELU_DEBUG
    // Test range of partition ids
    int incorrectRank = -1;
    for (int i = 0; i < decompEntries.size(); i++)
      if (decompEntries[i] >= numProcs || decompEntries[i] < 0) {
        incorrectRank = myRank;
        break;
      }

    int incorrectGlobalRank = -1;
    maxAll(comm, incorrectRank, incorrectGlobalRank);
    TEUCHOS_TEST_FOR_EXCEPTION(incorrectGlobalRank >- 1, Exceptions::RuntimeError, "pid " + toString(incorrectGlobalRank) + " encountered a partition number is that out-of-range");
#endif

    Array<GO> myGIDs;
    myGIDs.reserve(decomposition->getLocalLength());

    // Step 0: Construct mapping
    //    part number -> GIDs I own which belong to this part
    // NOTE: my own part GIDs are not part of the map
    typedef std::map<GO, Array<GO> > map_type;
    map_type sendMap;
    for (LO i = 0; i < decompEntries.size(); i++) {
      GO id  = decompEntries[i];
      GO GID = rowMap->getGlobalElement(i);

      if (id == myRank)
        myGIDs     .push_back(GID);
      else
        sendMap[id].push_back(GID);
    }

    if (keepProc0 && !remapPartitions) {
      // Figuring out how to keep processor 0 is easily and cheaply done in DeterminePartitionPlacement.
      // Here, we are in situation when DeterminePartitionPlacement was not called, but the user still
      // asked us to keep processor 0. Figuring that out is going to be slightly more difficult.

      // First, lets try to see if processor 0 gets any data. If it does, no need to do anything
      // For that, lets calculate the smalles part id that is valid.
      GO oldPartId, minPartId = Teuchos::OrdinalTraits<GO>::max();
      if (myGIDs.size())  minPartId = std::min(minPartId, Teuchos::as<GO>(myRank));
      if (sendMap.size()) minPartId = std::min(minPartId, sendMap.begin()->first);
      minAll(comm, minPartId, oldPartId);

      if (oldPartId == 0) {
        // Somebody owns a part with id 0. That means the processor 0 gets some data, even if it does
        // not have any originally. Our work is done.
        GetOStream(Statistics0, 0) << "No remapping is necessary despite that \"alwaysKeepProc0\" option is on,"
            " as processor 0 already receives some data" << std::endl;

      } else if (oldPartId == Teuchos::OrdinalTraits<GO>::max()) {
        // This is weird: nobody have any data. Nothing can be done.

      } else {
        // No partition with id 0, that means processor 0 gets no data. We have to do some extra legwork.
        // Specifically, we want to select a part such that the process owning the part has very small
        // fraction of the part.
        // NOTE: one could also trying minimizing the number of owned GIDs of that part, but assuming
        // good load balancing these metrics are the same.

        // Here is a neat trick: we can send minimizing information along with partition id but using
        // a single double. We use first numFracDigits digits of mantissa for actual fraction, and
        // numProcDigits digits after for storing the part id
        // NOTE: we need 10^{numAllDigits} to be smaller than INT_MAX
        const int    numFracDigits = 2,               numProcDigits = 7,               numAllDigits = numFracDigits + numProcDigits;
        const double powF = pow(10.0, numFracDigits), powP = pow(10.0, numProcDigits), powD = pow(10.0, numAllDigits);
        TEUCHOS_TEST_FOR_EXCEPTION(numProcs > powP, Exceptions::RuntimeError, "Time to update the constant!");

        double defaultFrac = 1.1, frac = defaultFrac, fracMin;
        if (myGIDs.size()) {
          frac = Teuchos::as<double>(myGIDs.size())/decompEntries.size();
        } else {
          // Some of the processors may have myGIDs size equal to zero. There are two way one could get that:
          //   1) Somebody sends pieces of part id = this processor id to it
          //   2) There is no part id corresponding to this processor id
          // Differentiatin between these two would require a lot more communication. Therefore, we exclude
          // all parts with no local GIDs from consideration. It results in suboptimal algorithm, but with
          // no extra communication
        }
        frac = (floor(frac*powF))/powF;                 // truncate the fraction to first numFracDigits
        frac = (floor(frac*powD) + myRank)/powD;        // store part id

        minAll(comm, frac, fracMin);

        if (fracMin < defaultFrac) {
          // Somebody sent some useful informtaion
          oldPartId = Teuchos::as<int>(fracMin*powD) % Teuchos::as<int>(powP); // decode

        } else {
          // Something weird is going on. This probably means that everybody does not keep any of its data
        }

        GetOStream(Statistics0, 0) << "Remapping part " << oldPartId << " to processor 0 as \"alwaysKeepProc0\" option is on" << std::endl;

        // Swap partitions
        // If a processor has a part of partition with id = oldPartId, that means that it sends data to it, unless
        // its rank is also oldPartId, in which case it some data is stored in myGIDs.
        if (myRank != 0 && myRank != oldPartId && sendMap.count(oldPartId)) {
          // We know for sure that there is no partition with id = 0 (there was a test for that). So we create one,
          // and swap the data with existing one.
          sendMap[0].swap(sendMap[oldPartId]);
          sendMap.erase(oldPartId);

        } else if (myRank == oldPartId && myGIDs.size()) {
          // We know for sure that there is no partition with id = 0 (there was a test for that). As all our data
          // belongs to processor 0 now, we move the data from myGIDs to the send array.
          sendMap[0].swap(myGIDs);

        } else if (myRank == 0 && sendMap.count(oldPartId)) {
          // We have some data that we send to oldPartId processor in the original distribution. We own that data now,
          // so we merge it with myGIDs array
          int offset = myGIDs.size(), len = sendMap[oldPartId].size();
          myGIDs.resize(offset + len);
          memcpy(myGIDs.getRawPtr() + offset, sendMap[oldPartId].getRawPtr(), len*sizeof(GO));
          sendMap.erase(oldPartId);
        }
      }
    }
    decompEntries = Teuchos::null;

    if (IsPrint(Statistics2)) {
      size_t numLocalKept = myGIDs.size(), numGlobalKept, numGlobalRows = A->getGlobalNumRows();
      sumAll(comm, numLocalKept, numGlobalKept);
      GetOStream(Statistics2,0) << "Unmoved rows: " << numGlobalKept << " / " << numGlobalRows << " (" << 100*Teuchos::as<double>(numGlobalKept)/numGlobalRows << "%)" << std::endl;
    }

    int numSend = sendMap.size(), numRecv;

    // Arrayify map keys
    Array<GO> myParts(numSend), myPart(1);
    int cnt = 0;
    myPart[0] = myRank;
    for (typename map_type::const_iterator it = sendMap.begin(); it != sendMap.end(); it++)
      myParts[cnt++] = it->first;

    // Step 1: Find out how many processors send me data
    RCP<Map>    partsIHave  = MapFactory   ::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), myParts(), indexBase, comm);
    RCP<Map>    partsIOwn   = MapFactory   ::Build(lib,                                                 numProcs,  myPart(), indexBase, comm);
    RCP<Export> partsExport = ExportFactory::Build(partsIHave, partsIOwn);

    RCP<GOVector> partsISend    = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(partsIHave);
    RCP<GOVector> numPartsIRecv = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(partsIOwn);
    if (numSend) {
      ArrayRCP<GO> partsISendData = partsISend->getDataNonConst(0);
      for (int i = 0; i < numSend; i++)
        partsISendData[i] = 1;
    }
    (numPartsIRecv->getDataNonConst(0))[0] = 0;

    numPartsIRecv->doExport(*partsISend, *partsExport, Xpetra::ADD);
    numRecv = (numPartsIRecv->getData(0))[0];

    // Step 2: Get my GIDs from everybody else
    MPI_Datatype MpiType = MpiTypeTraits<GO>::getType();
    int msgTag = 12345;  // TODO: use Comm::dup for all internal messaging

    // Post sends
    Array<MPI_Request> sendReqs(numSend);
    cnt = 0;
    for (typename map_type::iterator it = sendMap.begin(); it != sendMap.end(); it++)
      MPI_Isend(static_cast<void*>(it->second.getRawPtr()), it->second.size(), MpiType, Teuchos::as<GO>(it->first), msgTag, *rawMpiComm, &sendReqs[cnt++]);

    // Do waits
    map_type recvMap;
    size_t totalGIDs = myGIDs.size();
    for (int i = 0; i < numRecv; i++) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, msgTag, *rawMpiComm, &status);

      // Get rank and number of elements from status
      int fromRank = status.MPI_SOURCE, count;
      MPI_Get_count(&status, MpiType, &count);

      recvMap[fromRank].resize(count);
      MPI_Recv(static_cast<void*>(recvMap[fromRank].getRawPtr()), count, MpiType, fromRank, msgTag, *rawMpiComm, &status);

      totalGIDs += count;
    }

    // Merge GIDs
    myGIDs.reserve(totalGIDs);
    for (typename map_type::const_iterator it = recvMap.begin(); it != recvMap.end(); it++) {
      int offset = myGIDs.size(), len = it->second.size();
      if (len) {
        myGIDs.resize(offset + len);
        memcpy(myGIDs.getRawPtr() + offset, it->second.getRawPtr(), len*sizeof(GO));
      }
    }

    // Step 3: Construct importer
    RCP<Map>          newRowMap      = MapFactory   ::Build(lib, rowMap->getGlobalNumElements(), myGIDs(), indexBase, origComm);
    RCP<const Import> rowMapImporter = ImportFactory::Build(rowMap, newRowMap);

    Set(currentLevel, "Importer", rowMapImporter);

    // ======================================================================================================
    // Print some data
    // ======================================================================================================
    if (IsPrint(Statistics2)) {
      // Print the grid of processors
      GetOStream(Statistics2, 0) << "Partition distribution over cores (ownership is indicated by '+')" << std::endl;

      char amActive = (myGIDs.size() ? 1 : 0);
      std::vector<char> areActive(numProcs, 0);
      MPI_Gather(&amActive, 1, MPI_CHAR, &areActive[0], 1, MPI_CHAR, 0, *rawMpiComm);

      int rowWidth = Teuchos::as<int>(ceil(sqrt(numProcs)));
      for (int proc = 0; proc < numProcs; proc += rowWidth) {
        for (int j = 0; j < rowWidth; j++)
          if (proc + j < numProcs)
            GetOStream(Statistics2,0) << (areActive[proc + j] ? "+" : ".");
          else
          GetOStream(Statistics2,0) << " ";

        GetOStream(Statistics2,0) << "      " << proc << ":" << std::min(proc + rowWidth, numProcs) - 1 << std::endl;
      }
    }

  } // Build

  //----------------------------------------------------------------------
  template<typename T, typename W>
  struct Triplet {
    T    i, j;
    W    v;
  };
  template<typename T, typename W>
  static bool compareTriplets(const Triplet<T,W>& a, const Triplet<T,W>& b) {
    return (a.v > b.v); // descending order
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeterminePartitionPlacement(const Matrix& A, GOVector& decomposition,
                                                                                                               GO numPartitions, bool keepProc0) const {
    RCP<const Map>         rowMap    = A.getRowMap();
    GO                     indexBase = rowMap->getIndexBase();
    Xpetra::UnderlyingLib  lib       = rowMap->lib();

    RCP<const Teuchos::Comm<int> > comm = rowMap->getComm()->duplicate();
    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();


    // maxLocal is a constant which determins the number of largest edges which are being exchanged
    // The idea is that we do not want to construct the full bipartite graph, but simply a subset of
    // it, which requires less communication. By selecting largest local edges we hope to achieve
    // similar results but at a lower cost.
    const int maxLocal = 4;
    const int dataSize = 2*maxLocal;

    ArrayRCP<GO> decompEntries;
    if (decomposition.getLocalLength() > 0)
      decompEntries = decomposition.getDataNonConst(0);

    // Step 1: Sort local edges by weight
    // Each edge of a bipartite graph corresponds to a triplet (i, j, v) where
    //   i: processor id that has some piece of part with part_id = j
    //   j: part id
    //   v: weight of the edge
    // We set edge weights to be the total number of nonzeros in rows on this processor which
    // correspond to this part_id. The idea is that when we redistribute matrix, this weight
    // is a good approximation of the amount of data to move.
    // We use two maps, original which maps a partition id of an edge to the corresponding weight,
    // and a reverse one, which is necessary to sort by edges.
    std::map<GO,GO> lEdges;
    for (LO i = 0; i < decompEntries.size(); i++)
      lEdges[decompEntries[i]] += A.getNumEntriesInLocalRow(i);

    // Reverse map, so that edges are sorted by weight.
    // This results in multimap, as we may have edges with the same weight
    std::multimap<GO,GO> revlEdges;
    for (typename std::map<GO,GO>::const_iterator it = lEdges.begin(); it != lEdges.end(); it++)
      revlEdges.insert(std::make_pair(it->second, it->first));

    // Both lData and gData are arrays of data which we communicate. The data is stored
    // in pairs, so that data[2*i+0] is the part index, and data[2*i+1] is the corresponding edge weight.
    // We do not store processor id in data, as we can compute that by looking on the offset in the gData.
    Array<GO> lData(dataSize, -1), gData(numProcs * dataSize);
    int numEdges = 0;
    for (typename std::map<GO,GO>::reverse_iterator rit = revlEdges.rbegin(); rit != revlEdges.rend() && numEdges < maxLocal; rit++) {
      lData[2*numEdges+0] = rit->second; // part id
      lData[2*numEdges+1] = rit->first;  // edge weight
      numEdges++;
    }

    // Step 2: Gather most edges
    // Each processors contributes maxLocal edges by providing maxLocal pairs <part id, weight>, which is of size dataSize
    MPI_Datatype MpiType = MpiTypeTraits<GO>::getType();
    MPI_Allgather(static_cast<void*>(lData.getRawPtr()), dataSize, MpiType, static_cast<void*>(gData.getRawPtr()), dataSize, MpiType, *rawMpiComm);

    // Step 3: Construct mapping

    // Construct the set of triplets
    std::vector<Triplet<int,int> > gEdges(numProcs * maxLocal);
    size_t k = 0;
    for (LO i = 0; i < gData.size(); i += 2) {
      GO part   = gData[i+0];
      GO weight = gData[i+1];
      if (part != -1) {                     // skip nonexistent edges
        gEdges[k].i = i/dataSize;           // determine the processor by its offset (since every processor sends the same amount)
        gEdges[k].j = part;
        gEdges[k].v = weight;
        k++;
      }
    }
    gEdges.resize(k);

    // Sort edges by weight
    // NOTE: compareTriplets is actually a reverse sort, so the edges weight is in decreasing order
    std::sort(gEdges.begin(), gEdges.end(), compareTriplets<int,int>);

    // Do matching
    std::map<int,int> match;
    std::vector<char> matchedRanks(numProcs, 0), matchedParts(numProcs, 0);
    int lastMatchedPart = -1, numMatched = 0;
    for (typename std::vector<Triplet<int,int> >::const_iterator it = gEdges.begin(); it != gEdges.end(); it++) {
      GO rank = it->i;
      GO part = it->j;
      if (matchedRanks[rank] == 0 && matchedParts[part] == 0) {
        matchedRanks[rank] = 1;
        matchedParts[part] = 1;
        match[part] = rank;

        lastMatchedPart = part;
        numMatched++;
      }
    }
    GetOStream(Statistics0, 0) << "Number of unassigned paritions before cleanup stage: " << (numPartitions - numMatched) << " / " << numPartitions << std::endl;

    // Step 4 [optional]: Keep processor 0
    if (keepProc0)
      if (matchedRanks[0] == 0) {
        // Reassign partition to processor 0
        // The hope is that partition which we mapped last has few elements in it
        GetOStream(Statistics0, 0) << "Remapping part " << lastMatchedPart << " to processor 0 as \"alwaysKeepProc0\" option is on" << std::endl;
        matchedRanks[match[lastMatchedPart]] = 0;       // unassign processor which was matched to lastMatchedPart part
        matchedRanks[0] = 1;                            // assign processor 0
        match[lastMatchedPart] = 0;                     // match part to processor 0
      }

    // Step 5: Assign unassigned partitions
    // We do that through random matching for remaining partitions. Not all part numbers are valid, but valid parts are a subset of [0, numProcs).
    // The reason it is done this way is that we don't need any extra communication, as we don't need to know which parts are valid.
    for (int part = 0, matcher = 0; part < numProcs; part++)
      if (match.count(part) == 0) {
        // Find first non-matched rank
        while (matchedRanks[matcher])
          matcher++;

        match[part] = matcher++;
      }

    // Step 6: Permute entries in the decomposition vector
    for (LO i = 0; i < decompEntries.size(); i++)
      decompEntries[i] = match[decompEntries[i]];
  }

} // namespace MueLu

#endif //ifdef HAVE_MPI

#endif // MUELU_REPARTITIONFACTORY_DEF_HPP
