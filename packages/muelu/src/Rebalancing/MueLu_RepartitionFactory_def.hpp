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

#include "MueLu_Utilities.hpp" // TMP JG NOTE: only for maxAll, so no _fwd in _decl

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<int>        ("startLevel",                   1, "First level at which repartitioning can possibly occur. Repartitioning at finer levels is suppressed");
    validParamList->set<LO>         ("minRowsPerProcessor",       1000, "Minimum number of rows over all processes. If any process falls below this, repartitioning is initiated");
    validParamList->set<double>     ("nonzeroImbalance",           1.2, "Imbalance threshold, below which repartitioning is initiated. Imbalance is measured by "
                                                                        "ratio of maximum nonzeros over all processes to minimum number of nonzeros over all processes");
    validParamList->set<bool>       ("fixedOrder",                true, "Use sorting of recv PIDs to force reproducibility");
    // FIXME: Unused; LO instead of GO?
    // validParamList->set<GO>         ("minNnzPerProcessor",          -1, "Minimum number of nonzeros over all processes. If any process falls below this, repartitioning is initiated.");
    validParamList->set<std::string>("adjustNumPartitions",     "none", "Algorithm for adjusting number of partitions (none|2k)");

    {
      std::stringstream docDiffusiveHeuristic;
      docDiffusiveHeuristic << "0:   put on procs 0..N" << std::endl
                            << "1:   use diffusive heuristic" << std::endl
                            << "K>1: if #partitions is > K, put on procs 0..N, otherwise use diffusive heuristic" << std::endl
                            << "-1:  put on procs 0..N the first time only, then use diffusive in remaining rounds" << std::endl;

      validParamList->set<LO>("diffusiveHeuristic",     0, docDiffusiveHeuristic.str());
    }
    validParamList->set<LO>("maxNumArbitrationRounds",     10, "Maximum number of arbitration rounds for diffusive heuristic");

    validParamList->set< RCP<const FactoryBase> >("A",                   Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Partition",           Teuchos::null, "Factory of the partition");

    // One can specify a number of partition on the level class to override what is computed internally by this class.
    // By default, the generation factory is 'this', not Teuchos::null, which means that by default, neither user defined data nor
    // the factory manager are used. This is unusual.
    // To manually set a "number of partition" entry in the level, it has to either be associated with the generating factory 'this' or the generating factory of this class has to be set to Teuchos::null:
    // Ex: level.Set("numbers of partitions"); myRepartitionFact.set< RCP<FactoryBase> >("number of partitions", Teuchos::null);
    //  or level.Set("numbers of partitions", myRepartitionFact); and approrpiate requests
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

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    const Teuchos::ParameterList & pL = GetParameterList();
    // Access parameters here to make sure that we set the parameter entry flag to "used" even in case of short-circuit evaluation.
    // JG TODO: I don't really know if we want to do this.
    const int    startLevel          = pL.get<int>   ("startLevel");
    const LO     minRowsPerProcessor = pL.get<LO>    ("minRowsPerProcessor");
    const double nonzeroImbalance    = pL.get<double>("nonzeroImbalance");

    //TODO: We only need a CrsGraph. This class does not have to be templated on Scalar types.
    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

    { // scoping

      // ======================================================================================================
      // Determine whether partitioning is needed.
      // ======================================================================================================
      // TODO: Result of the test can be stored in the level to avoid doing this several time.

      // Repartitionning iff !(Test1 || Test2) && (Test3 || Test4)
      // This is implemented with a short-circuit evaluation

      // Test1: skip partitioning if level is too big
      std::ostringstream msg1; msg1 << std::endl << "    current level = " << currentLevel.GetLevelID() << ", first level where repartitioning can happen is " << startLevel << ".";
      if (currentLevel.GetLevelID() < startLevel) {
        GetOStream(Warnings0, 0) << "No repartitioning necessary:" + msg1.str() << std::endl;
        Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
        return;
      }

      // Test 2: check whether A is spread over more than one process.
      // Note: using type 'int' as it is enough (comm->getSize() is an int). Test can also be done with 'bool' but numActiveProcesses is printed later
      // TODO: this global communication can be avoided if we store the information with the matrix (it is known when matrix is created)
      // TODO: further improvements could be achieved when we use subcommunicator for the active set. Then we only need to check its size
      int numActiveProcesses = 0;
      sumAll(comm, (int)((A->getNodeNumRows() > 0) ? 1 : 0), numActiveProcesses);
      std::ostringstream msg2; msg2 << std::endl << "    # processes with rows = " << numActiveProcesses;
      if (numActiveProcesses == 1) {
        GetOStream(Warnings0, 0) << "No repartitioning necessary:" + msg2.str() << std::endl;
        Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
        return;
      }

      bool doRepartition = false;
      std::ostringstream msg3, msg4, *msgDoRepartition=NULL; // msgDoRepartition == msg3 or msg4 depending on Test3 and Test4 evalution

      // Test3: check whether any node has too few rows
      // Note: (!Test2) ensures that repartitionning is not done when only 1 proc and globalNumRow < minRowsPerProcessor
      // TODO: change it to repartition if only the number of processors with numRows < threshold is larger than some percentage
      // of the total number. This way, we won't repartition if 2 out of 1000 processors don't have enough elements. I'm thinking
      // maybe 20% threshold.
      if (minRowsPerProcessor > 0) { // skip test if criteria not used
        LO     minNumRows;
        size_t numMyRows = A->getNodeNumRows();
        LO LOMAX = Teuchos::OrdinalTraits<LO>::max(); // processors without rows do not participate to minAll()
        minAll(comm, (LO)((numMyRows > 0) ? numMyRows : LOMAX), minNumRows);
        TEUCHOS_TEST_FOR_EXCEPTION(minNumRows >= LOMAX, Exceptions::RuntimeError, "internal error");
        msg3 << std::endl << "    min # rows per proc = "   << minNumRows << ", min allowable = " << minRowsPerProcessor;
        if (minNumRows < minRowsPerProcessor) {
          doRepartition = true;
          msgDoRepartition = &msg3;
        }
      }

      // Test4: check whether the number of nonzeros per process is imbalanced
      if (doRepartition == false) { // skip test if we already now that we will do repartitioning
        size_t numMyNnz  = A->getNodeNumEntries();
        GO minNnz, maxNnz;
        maxAll(comm, (GO)numMyNnz, maxNnz);
        minAll(comm, (GO)((numMyNnz > 0) ? numMyNnz : maxNnz), minNnz); // min nnz over all proc (disallow any processors with 0 nnz)
        double imbalance = ((double) maxNnz) / minNnz;
        msg4 << std::endl << "    nonzero imbalance = " << imbalance  << ", max allowable = " << nonzeroImbalance;
        if (imbalance > nonzeroImbalance) {
          doRepartition = true;
          msgDoRepartition = &msg4;
        }
      }

      if (!doRepartition) {
        GetOStream(Warnings0, 0) << "No repartitioning necessary:" + /* msg1.str() + msg2.str() +*/ msg3.str() + msg4.str() << std::endl;
        Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
        return;
      }

      // print only conditions that triggered the repartitioning
      GetOStream(Statistics0, 0) << "Repartitioning necessary:" << msg1.str() << msg2.str() << msgDoRepartition->str() << std::endl;

    } // scoping

    // FIXME Quick way to figure out how many partitions there should be. (Same as what's done in ML.)
    // FIXME Should take into account nnz, maybe?  Perhaps only when user is using min #nnz per row threshold.
    GO numPartitions;
    if (IsAvailable(currentLevel, "number of partitions")) {
      numPartitions = Get<GO>(currentLevel, "number of partitions");

    } else {
      if ((GO)A->getGlobalNumRows() < minRowsPerProcessor) {
        // System is too small, migrate it to a single processor
        numPartitions = 1;

      } else {
        // Make sure that each processor has approximately minRowsPerProcessor
        numPartitions = A->getGlobalNumRows() / minRowsPerProcessor;
      }
      if (numPartitions > comm->getSize())
        numPartitions = comm->getSize();

      std::string adjustment = pL.get<std::string>("adjustNumPartitions");
      if (adjustment == "2k") {
        GetOStream(Statistics0, 0) << "Number of partitions to use = " << numPartitions << std::endl;

        int i2 = Teuchos::as<int>(floor(log(numPartitions)/log(2)));
        numPartitions = Teuchos::as<int>(std::pow(2.,i2));

        GetOStream(Runtime0,0) << "Adjusting number of partitions using \"2k\" algorithm to " << numPartitions << std::endl;
      }

      Set(currentLevel, "number of partitions", numPartitions);
    }
    GetOStream(Statistics0, 0) << "Number of partitions to use = " << numPartitions << std::endl;

    // ======================================================================================================
    // Determine the global size of each partition.
    // ======================================================================================================
    // Length of vector "decomposition" is local number of DOFs.  Its entries are partition numbers each DOF belongs to.

    RCP<GOVector> decomposition;
    if (numPartitions == 1) {
      // Running on one processor, so decomposition is the trivial one, all zeros.
      // => skip call to Zoltan_Interface (There is also a short cut in Zoltan_Interface so this is mainly to avoid extra output messages)
      // TODO: can we skip more work (ie: building the hashtable, etc.)?
      GetOStream(Warnings0, 0) << "Only one partition: Skip call to the repartitioner." << std::endl;
      decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), true);
    } else {
      decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");

      // Zoltan2 changes the number of partitions. There is no good mechanism to propagate that new number
      // to this factory, but we can do that by finding out the max number of partition across all processors
      GO maxPartLocal = 0;
      if (decomposition->getLocalLength()) {
        // NOTE: this is a stupid check. We would not be here if we didn't have any data.
        // But one of the unit tests (Repartition_Build) constructs a stupid map in which processor
        // 2 does not have any data. I'll add this check here.
        maxPartLocal = *std::max_element(decomposition->getData(0).begin(), decomposition->getData(0).end());
      }
      maxAll(decomposition->getMap()->getComm(), maxPartLocal, numPartitions);
      numPartitions++;

      Set(currentLevel, "number of partitions", numPartitions);

      if (decomposition == Teuchos::null) {
        GetOStream(Warnings0, 0) << "No repartitioning necessary: partitions were left unchanged by the repartitioner" << std::endl;
        Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
        return;
      }
    }

    // Use a hashtable to record how many local rows belong to each partition.
    RCP<Teuchos::Hashtable<GO, GO> > hashTable;
    hashTable = rcp(new Teuchos::Hashtable<GO, GO>(numPartitions + numPartitions/2));
    Teuchos::Hashtable<GO,GO>& htref = *hashTable;
    ArrayRCP<const GO> decompEntries;
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    bool flag=false;
    for (int i=0; i<decompEntries.size(); ++i) {
      if (decompEntries[i] >= numPartitions) flag = true;
      if (htref.containsKey(decompEntries[i])) {
        GO count = htref.get(decompEntries[i]);
        ++count;
        htref.put(decompEntries[i], count);
      } else {
        htref.put(decompEntries[i], 1);
      }
    }
    int problemPid;
    int mypid = comm->getRank();
    maxAll(comm, (flag ? mypid : -1), problemPid);
    std::ostringstream buf; buf << problemPid;
    TEUCHOS_TEST_FOR_EXCEPTION(problemPid>-1, Exceptions::RuntimeError, "pid " + buf.str() + " encountered a partition number is that out-of-range");
    decompEntries = Teuchos::null;

    Teuchos::Array<GO> allPartitionsIContributeTo;
    Teuchos::Array<GO> allLocalPartSize;
    htref.arrayify(allPartitionsIContributeTo, allLocalPartSize);

    GO indexBase = decomposition->getMap()->getIndexBase();

    // Source map is overlapping.  GIDs owned by this pid are the partition numbers found above.
    RCP<Map> sourceMap = Xpetra::MapFactory<LO, GO, NO>::Build(decomposition->getMap()->lib(),
                                           Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                           allPartitionsIContributeTo(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

    // Store # of local DOFs in each partition in a vector based on above map.
    RCP<GOVector> localPartSizeVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(sourceMap, false);
    ArrayRCP<GO> data;
    if (localPartSizeVec->getLocalLength() > 0)
      data = localPartSizeVec->getDataNonConst(0);
    for (int i=0; i<htref.size(); ++i)
      data[i] = allLocalPartSize[i];
    data = Teuchos::null;

    // Target map is nonoverlapping.  Pid k has GID N if and only if k owns partition N.
    GO myPartitionNumber;
    Array<int> partitionOwners;
    partitionOwners.reserve(numPartitions);
    DeterminePartitionPlacement(currentLevel, myPartitionNumber, partitionOwners);

/*
    // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
    std::cout << "pid " << mypid << " here" << std::endl;
    sleep(1); comm->barrier();
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << " owns partition " << myPartitionNumber << std::endl;
        for (int j=0; j<partitionOwners.size(); ++j)
          std::cout << "     partition " << j << " owned by pid " << partitionOwners[j] << std::endl;
      }
      comm->barrier();
    }
    sleep(1); comm->barrier();
    // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
*/

    GO numDofsThatStayWithMe=0;
    Teuchos::Array<GO> partitionsIContributeTo;
    partitionsIContributeTo.reserve(allPartitionsIContributeTo.size());
    Teuchos::Array<GO> localPartSize;
    localPartSize.reserve(allPartitionsIContributeTo.size());
    for (int i=0; i<allPartitionsIContributeTo.size(); ++i)
    {
      if (allPartitionsIContributeTo[i] != myPartitionNumber) {
        partitionsIContributeTo.push_back(allPartitionsIContributeTo[i]);
        localPartSize.push_back(allLocalPartSize[i]);
      }
      else
        numDofsThatStayWithMe = allLocalPartSize[i];
    }

    // Note: "numPartitionsISendTo" does not include this PID
    GO numPartitionsISendTo = htref.size();
    // I'm a partition owner, so don't count my own.
    if (myPartitionNumber >= 0 && numPartitionsISendTo>0 && numDofsThatStayWithMe>0) numPartitionsISendTo--;
    assert(numPartitionsISendTo == partitionsIContributeTo.size());

    Array<int> partitionOwnersISendTo;
    partitionOwnersISendTo.reserve(allPartitionsIContributeTo.size());
    for (int i=0; i<partitionsIContributeTo.size(); ++i) {
      partitionOwnersISendTo.push_back(partitionOwners[partitionsIContributeTo[i]]);
    }

    Array<GO> localMapElement;
    if (myPartitionNumber >= 0)
      localMapElement.push_back(myPartitionNumber);
    RCP<Map> targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           numPartitions,
                                           localMapElement(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

    RCP<const Export> exporter = ExportFactory::Build( sourceMap, targetMap);

    // If this pid owns a partition, globalPartSizeVec has one local entry that is the global size of said partition.
    // If this pid doesn't own a partition, globalPartSizeVec is locally empty.
    RCP<GOVector> globalPartSizeVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(targetMap);
    globalPartSizeVec->doExport(*localPartSizeVec, *exporter, Xpetra::ADD);
    int myPartitionSize = 0;
    ArrayRCP<const GO> constData;
    if (globalPartSizeVec->getLocalLength() > 0) {
      constData = globalPartSizeVec->getData(0);
      myPartitionSize = constData[0];
    }
    constData = Teuchos::null;

    // ======================================================================================================
    // Calculate how many PIDs (other than myself) contribute to my partition.
    // ======================================================================================================
    RCP<GOVector> howManyPidsSendToThisPartitionVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(targetMap);
    RCP<GOVector> partitionsISendTo = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(sourceMap, false);
    if (partitionsISendTo->getLocalLength() > 0)
      data = partitionsISendTo->getDataNonConst(0);
    const Map& sourceMapRef = *sourceMap;
    for (int i=0; i<data.size(); ++i) {
      // don't count myself as someone I send to... (sourceMap is based on allPartitionsIContributeTo)
      if (sourceMapRef.getGlobalElement(i) != myPartitionNumber)
        data[i] = 1;
      else
        data[i] = 0;
    }
    data = Teuchos::null;

    // Note: "howManyPidsSendToThisPartition" does not include this PID
    int howManyPidsSendToThisPartition = 0;
    howManyPidsSendToThisPartitionVec->doExport(*partitionsISendTo, *exporter, Xpetra::ADD);
    if (howManyPidsSendToThisPartitionVec->getLocalLength() > 0) {
      constData = howManyPidsSendToThisPartitionVec->getDataNonConst(0);
      howManyPidsSendToThisPartition = constData[0];
    }
    constData = Teuchos::null;

    // ======================================================================================================
    // Calculate which PIDs contribute to my partition, and how much each contributes,
    // with ireceive/send/wait cycle.
    // FIXME Jan.12.2012 Teuchos::Comm methods ireceive and wait don't work (bugs 5483 and 5484), so for
    // now use the raw MPI methods.
    // ======================================================================================================

    Array<GO> pidsIReceiveFrom(howManyPidsSendToThisPartition);
    Array<GO> numDofsIReceiveFromOnePid(howManyPidsSendToThisPartition);
    for (int j=0; j<numDofsIReceiveFromOnePid.size(); ++j)
      numDofsIReceiveFromOnePid[j] = -99;

    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    if (tmpic == Teuchos::null)
      throw(Exceptions::RuntimeError("Cannot cast base Teuchos::Comm to Teuchos::MpiComm object."));
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

    // First post non-blocking receives.
    Array<MPI_Request> requests(howManyPidsSendToThisPartition);
    for (int i=0; i<howManyPidsSendToThisPartition; ++i) {
      MPI_Irecv((void*)&(numDofsIReceiveFromOnePid[i]), sizeof(GO), MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, *rawMpiComm, &requests[i]);
    }

    // Next post sends.
    for (int i=0; i< partitionsIContributeTo.size(); ++i) {
      comm->send(sizeof(GO), (char*)&localPartSize[i], partitionOwnersISendTo[i]);
    }

    // Finally do waits.
    Array<MPI_Status> status(howManyPidsSendToThisPartition);
    for (int i=0; i<howManyPidsSendToThisPartition; ++i)
      MPI_Wait(&requests[i], &status[i]);

    for (int i=0; i<pidsIReceiveFrom.size(); ++i)
      pidsIReceiveFrom[i] = status[i].MPI_SOURCE;

    if (pL.get<bool>("fixedOrder")) {
      // Sort pids for reproducibility (and numDofsIReceiveFromOnePid for consistency)
      GetOStream(Runtime0,0) << "Sorting recv PIDs for results reproducibility" << std::endl;

      // For simplicity, do insertion sort of pidsIReceiveFrom and use that permutation for numDofsIReceiveFrom
      // NOTE: if insertion sort becomes slow for some use cases, we may replace it with the quicksort Tpetra::sort2
      for (int i = 1; i < howManyPidsSendToThisPartition; i++) {
        GO d = pidsIReceiveFrom[i], v = numDofsIReceiveFromOnePid[i];
        int j;
        for (j = i-1; j >= 0 && pidsIReceiveFrom[j] > d; j--) {
          pidsIReceiveFrom         [j+1] = pidsIReceiveFrom[j];
          numDofsIReceiveFromOnePid[j+1] = numDofsIReceiveFromOnePid[j];
        }
        pidsIReceiveFrom         [j+1] = d;
        numDofsIReceiveFromOnePid[j+1] = v;
      }
    }

    comm->barrier();

    // =================================================================================================
    // Calculate partial offsets for permutation row map, via MPI_Scan based on global partition sizes.
    // Communicate those offsets back to respective PIDS of unpermuted matrix using ireceive/send/wait cycle.
    // =================================================================================================
    // Partition numbers are used as tags to ensure messages go to correct PIDs.
    // Use MPI_DOUBLE type to avoid any overflow problems.
    // partitionSizeOffset is the first GID in this partition.
    // gidOffsets is an array of GID offsets.  gidOffsets[i] is the first GID for the dofs received from PID partitionOwnersISendTo[i].
    // Note: Quantities "numPartitionsISendTo" and "howManyPidsSendToThisPartition" do not include me
    double partitionSizeOffset;
    double ttt = myPartitionSize;
    MPI_Scan(&ttt, &partitionSizeOffset, 1, MPI_DOUBLE, MPI_SUM, *rawMpiComm);
    partitionSizeOffset -= myPartitionSize;

    Array<double> gidOffsets;
    if (howManyPidsSendToThisPartition > 0) {
      gidOffsets.resize(howManyPidsSendToThisPartition);
      gidOffsets[0] = partitionSizeOffset;
    }
    for (int i=1; i<howManyPidsSendToThisPartition; ++i) {
      gidOffsets[i] = gidOffsets[i-1] + numDofsIReceiveFromOnePid[i-1];
    }

    requests.resize(numPartitionsISendTo);
    Array<double> gidOffsetsForPartitionsIContributeTo(numPartitionsISendTo);

    // Post receives on contributing PIDs.
    for (int i=0; i< numPartitionsISendTo; ++i) {
      int msgTag = partitionsIContributeTo[i];
      MPI_Irecv((void*)&(gidOffsetsForPartitionsIContributeTo[i]), 1, MPI_DOUBLE, partitionOwnersISendTo[i], msgTag, *rawMpiComm, &requests[i]);
    }

    // Do sends by partition owners.
    for (int i=0; i<howManyPidsSendToThisPartition; ++i) {
      int msgTag = myPartitionNumber;
      MPI_Send((void*)&gidOffsets[i] , 1 , MPI_DOUBLE , pidsIReceiveFrom[i], msgTag, *rawMpiComm);
    }

    // Do waits.
    status.resize(numPartitionsISendTo);
    for (int i=0; i<numPartitionsISendTo; ++i)
      MPI_Wait(&requests[i], &status[i]);

    comm->barrier();

    // =================================================================================================
    // Set up a synthetic GID scheme that is the same for both the original unpermuted system and the permuted system.
    // This scheme is *not* the final DOF numbering, but is just for the importer we need to transfer column IDs.
    // =================================================================================================

    // Synthetic GIDS for original unpermuted system.

    // store offsets for easy random access
    hashTable = rcp(new Teuchos::Hashtable<GO, GO>(partitionsIContributeTo.size() + partitionsIContributeTo.size()/2));
    Teuchos::Hashtable<GO,GO>& htref1 = *hashTable;
    for (int i=0; i<partitionsIContributeTo.size(); ++i) {
      htref1.put(partitionsIContributeTo[i], (GO)gidOffsetsForPartitionsIContributeTo[i]);
    }
    //store gid offset for those dofs that will remain with me
    if (myPartitionNumber > -1)
      htref1.put(myPartitionNumber, ((GO)partitionSizeOffset) + myPartitionSize - numDofsThatStayWithMe);
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    Array<GO> uniqueGIDsBeforePermute;
    uniqueGIDsBeforePermute.reserve(decompEntries.size());
    for (int i=0; i<decompEntries.size(); ++i) {
      GO gid = htref1.get(decompEntries[i]);
      uniqueGIDsBeforePermute.push_back(gid);
      gid++;
      htref1.put(decompEntries[i], gid);
    }
    decompEntries = Teuchos::null;

    // Synthetic GIDS for permuted system.

    Array<GO> uniqueGIDsAfterPermute;
    uniqueGIDsAfterPermute.reserve(myPartitionSize);
    for (int i=0; i<myPartitionSize; ++i) {
        uniqueGIDsAfterPermute.push_back((GO)partitionSizeOffset+i);
    }

    // =================================================================================================
    // Create and apply an importer to communicate row GIDs that will make up the permuted row map
    // =================================================================================================

    //TODO we should really supply the global size as another sanity check
    sourceMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    uniqueGIDsBeforePermute(),
                    indexBase,
                    comm);

    targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    uniqueGIDsAfterPermute(),
                    indexBase,
                    comm);

    RCP<const Import> importer = ImportFactory::Build( sourceMap, targetMap);

    RCP<GOVector> sourceVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(sourceMap);
    ArrayRCP<GO> vectorData;
    if (sourceVec->getLocalLength() > 0) {
      vectorData = sourceVec->getDataNonConst(0);
    }

    //load source vector with unpermuted GIDs.
    //RCP<const Map> originalRowMap = A->getRowMap();
    const Map& originalRowMap = *(A->getRowMap());

    //sanity check
    assert(vectorData.size() == (GO) originalRowMap.getNodeNumElements());


    for (LO i=0; i<vectorData.size(); ++i) {
      GO gid = originalRowMap.getGlobalElement(i);
      if (gid == (GO) Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid())
        throw(Exceptions::RuntimeError("Encountered an unexpected error with A's rowmap."));
      vectorData[i] = gid;
    }
    vectorData = Teuchos::null;

    RCP<GOVector> targetVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(targetMap);
    targetVec->doImport(*sourceVec, *importer, Xpetra::INSERT);

    // =================================================================================================
    // Create an importer between the original row map and the permuted row map
    // =================================================================================================
    vectorData = targetVec->getDataNonConst(0);
    RCP<Map> newRowMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    vectorData(),
                    indexBase,
                    comm);

    RCP<const Import> importerForRepartitioning = ImportFactory::Build(A->getRowMap(), newRowMap);
    Set(currentLevel, "Importer", importerForRepartitioning);

  } //Build

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeterminePartitionPlacement(Level & currentLevel, GO &myPartitionNumber,
  Array<int> &partitionOwners) const
  {
    FactoryMonitor m(*this, "DeterminePartitionPlacement", currentLevel);

    GO numPartitions = Get<GO>(currentLevel, "number of partitions");
    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
    int mypid = comm->getRank();

    const Teuchos::ParameterList & pL = GetParameterList();
    const LO diffusiveHeuristic = pL.get<LO>("diffusiveHeuristic");
    static bool forceDiffusive = false; // FIXME: we need a mecanism to reset this variable (if several Hierarchy::Setup() using same factories)

    // If not diffusive, compute (partitionOwners, myPartitionNumber) and return.
    if (diffusiveHeuristic != 1 && numPartitions > diffusiveHeuristic && !forceDiffusive) {

      // if diffusiveHeuristic == -1:  put on procs 0..N the first time only, then use diffusive in remaining rounds
      if (diffusiveHeuristic == -1) forceDiffusive = true;

      if (numPartitions==1)
        GetOStream(Runtime0, 0) << "Placing partitions on proc. 0." << std::endl;
      else
        GetOStream(Runtime0, 0) << "Placing partitions on proc. 0-" << numPartitions-1 << "." << std::endl;

      myPartitionNumber = -1;
      for (int i=0; i<(int)numPartitions; ++i) {
        partitionOwners.push_back(i);
        if (i==mypid) myPartitionNumber = i;
      }

      return;
    }

    GetOStream(Runtime0, 0) << "Using diffusive heuristic for partition placement." << std::endl;

    RCP<SubFactoryMonitor> m1 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: Setup", currentLevel));

//RCP<SubFactoryMonitor> m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: getting 'Partition'", currentLevel));
    RCP<GOVector> decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");
    // Figure out how many nnz there are per row.
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: figuring out nnz per row", currentLevel));
    RCP<GOVector> nnzPerRowVector = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO> nnzPerRow;
    if (nnzPerRowVector->getLocalLength() > 0)
      nnzPerRow = nnzPerRowVector->getDataNonConst(0);
    for (int i=0; i<nnzPerRow.size(); ++i)
      nnzPerRow[i] = A->getNumEntriesInLocalRow(i);

    // Use a hashtable to record how many nonzeros in the local matrix belong to each partition.
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: hashing", currentLevel));
    RCP<Teuchos::Hashtable<GO, GO> > hashTable;
    hashTable = rcp(new Teuchos::Hashtable<GO, GO>(numPartitions + numPartitions/2));
    Teuchos::Hashtable<GO, GO>& htref = *hashTable;

    //decomposition is a vector such that decomposition[i] is the partition that local row i belongs to
    //numPartitions is the global # of partitions
    ArrayRCP<const GO> decompEntries;
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    bool flag=false;
    assert(decompEntries.size() == nnzPerRow.size());
    for (int i=0; i<decompEntries.size(); ++i) {
      if (decompEntries[i] >= numPartitions) flag = true; //error checking (see exception test below)
      if (htref.containsKey(decompEntries[i])) {
        GO count = htref.get(decompEntries[i]);
        count += nnzPerRow[i];
        htref.put(decompEntries[i], count);
      } else {
        htref.put(decompEntries[i], nnzPerRow[i]);
      }
    }

//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: arrayify", currentLevel));
    //error checking FIXME skip the maxall and just throw already
    int problemPid;
    maxAll(comm, (flag ? mypid : -1), problemPid);
    std::ostringstream buf; buf << problemPid;
    TEUCHOS_TEST_FOR_EXCEPTION(problemPid>-1, Exceptions::RuntimeError, "pid " + buf.str() + " encountered a partition number is that out-of-range");
    decompEntries = Teuchos::null;

    Teuchos::Array<GO> allPartitionsIContributeTo;
    Teuchos::Array<GO> localNnzPerPartition;
    htref.arrayify(allPartitionsIContributeTo, localNnzPerPartition);
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build target map", currentLevel));
    //map in which all pids have all partition numbers as GIDs.
    //this map has duplicated GIDs, obviously
    //FIXME this next map ctor can be a real time hog in parallel
    /*
    Array<GO> allPartitions;
    allPartitions.reserve(numPartitions);
    for (int i=0; i<numPartitions; ++i) allPartitions.push_back(i);
    RCP<Map> targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           numPartitions*comm->getSize(),
                                           allPartitions(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);
    */
    // Locally, this map has only partition GIDs for which the PID has dofs for
    RCP<Map> targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                           allPartitionsIContributeTo(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build vectors", currentLevel));
    RCP<Xpetra::Vector<double, LO, GO, NO> > globalWeightVec = Xpetra::VectorFactory<double, LO, GO, NO>::Build(targetMap);  //TODO why does the compiler grumble about this when I omit template arguments?
    RCP<Xpetra::Vector<LO, LO, GO, NO> > procWinnerVec = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(targetMap);
    ArrayRCP<LO> procWinner;
    if (procWinnerVec->getLocalLength() > 0)
      procWinner = procWinnerVec->getDataNonConst(0);
    for (int i=0; i<procWinner.size(); ++i) procWinner[i] = -1;
    procWinner = Teuchos::null;
    RCP<Xpetra::Vector<SC, LO, GO, NO> > scalarProcWinnerVec = Xpetra::VectorFactory<SC, LO, GO, NO>::Build(targetMap);

//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build unique map", currentLevel));
    Array<GO> myPidArray;
    myPidArray.push_back(mypid);

    // intermediate map required by ArbitrateAndCommunicate
    RCP<Map> uniqueMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           comm->getSize(),
                                           myPidArray(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build comm helper", currentLevel));
    MueLu::CoupledAggregationCommHelper<LO, GO, NO, LMO> commHelper(uniqueMap, targetMap);
    myPartitionNumber = -1;
    int doArbitrate = 1;
//m3 = Teuchos::null;

    /*
       Use ArbitrateAndCommunicate to determine which process should own each partition.
       This may require multiple rounds because a single process i may be found to be the
       largest contributor for more than one partition (say P1, P2, ...Pn).  If that happens,
       process i is made owner of the partition Pk to which it contributes the most.  If
       two processes end up wanting to own the same partition, it's the job of A&C to break
       the tie.

       The contributions of process i to all other partitions are set to 0 (effectively meaning
       i can't be assigned another partition by A&C), and A&C is called again.

       continues until all partitions have been uniquely assigned.
    */

    m1 = Teuchos::null;

    int numRounds=0;
    RCP<SubFactoryMonitor> m2 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: Arbitration phase", currentLevel));

    const LO maxRounds = pL.get<LO>("maxNumArbitrationRounds");

    while (doArbitrate && numRounds < maxRounds)
    {
      ++numRounds;
      ArrayRCP<double> globalWeightVecData = globalWeightVec->getDataNonConst(0);

      //If this process doesn't yet own a partition, record all its nonzeros per partition as weights.
      //If it doesn't contribute to a partition, make the weight zero.
      if (myPartitionNumber == -1) {
        for (int i=0; i<globalWeightVecData.size(); ++i)
          globalWeightVecData[i] = 0.0;
        for (int i=0; i<allPartitionsIContributeTo.size(); ++i)
          globalWeightVecData[ targetMap->getLocalElement(allPartitionsIContributeTo[i]) ] = localNnzPerPartition[i];
      } else {
        //this process already owns a partition, so record only the #nonzeros in that partition as weights
        for (int i=0; i<globalWeightVecData.size(); ++i)
          globalWeightVecData[i] = 0.0;
        for (int i=0; i<allPartitionsIContributeTo.size(); ++i) {
          if (allPartitionsIContributeTo[i] == myPartitionNumber) {
            globalWeightVecData[ targetMap->getLocalElement(myPartitionNumber) ] = localNnzPerPartition[i];
            break;
          }
        }
      }

      globalWeightVecData = Teuchos::null;

      commHelper.ArbitrateAndCommunicate(*globalWeightVec, *procWinnerVec, NULL, true);

      procWinner = procWinnerVec->getDataNonConst(0);

      // If this process has tentatively been assigned more than one partition,
      // choose the largest as the partition that this process will own.
      // Due to how ArbitrateAndCommunicate works, do the following:
      //   + Mark as unowned all other smaller partitions that were tentatively assigned to me.
      //   + Mark as unowned all partitions that I do not own.
      GO partitionsIOwn=0;
      GO largestPartitionSize=0;

      int lastWon=-1;
      for (int i=0; i<procWinner.size(); ++i) {
        if (procWinner[i] == mypid) {
          GO partitionSize = htref.get(targetMap->getGlobalElement(i));
          if (partitionSize > largestPartitionSize) {
            myPartitionNumber = targetMap->getGlobalElement(i);
            largestPartitionSize = partitionSize;
            if (lastWon>-1)
              procWinner[lastWon] = -1;
            lastWon = i;
          } else
            procWinner[i] = -1;
          partitionsIOwn++;
        } else {
          procWinner[i] = -1;
        }
     }
      procWinner = Teuchos::null;

      // If this pid was tentatively assigned more than one partition, then we can arbitrate again
      // because there are unassigned partitions.
      int arbitrateAgain;
      if (partitionsIOwn > 1) arbitrateAgain=1;
      else                    arbitrateAgain=0;
      sumAll(comm, arbitrateAgain, doArbitrate);

#     ifdef HAVE_MUELU_DEBUG
      GetOStream(Statistics0, 0) << "round " << numRounds << " : ";
      GO ii=0, jj=0;
      if (partitionsIOwn == 1) ii=1;
      sumAll(comm, ii,jj);
      GetOStream(Statistics0, 0) << "#definitely assigned = " << jj;
      ii=0, jj=0;
      if (partitionsIOwn > 1) ii=1;
      sumAll(comm, ii, jj);
      GetOStream(Statistics0, 0) << ", #pids with > 1 partition =  " << jj;
      jj=0;
      maxAll(comm, partitionsIOwn, jj);
      GetOStream(Statistics0, 0) << ", max# partitions owned by single pid = " << jj;
      GO kk=0; ii=0;
      if (partitionsIOwn == jj) kk=1;
      sumAll(comm, kk, ii);
      if (ii>1)
        GetOStream(Statistics0, 0) << "(by " << ii << " procs)";
      if (ii==1) {
        ii=0; kk=0;
        if (partitionsIOwn == jj) kk=Teuchos::as<GO>(comm->getRank());
        sumAll(comm, kk, ii);
        GetOStream(Statistics0, 0) << "(by pid " << ii << ")";
      }
      comm->barrier();
      GetOStream(Statistics0, 0) << std::endl;
#     endif

    } //while (doArbitrate)
    m2 = Teuchos::null;

    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    if (tmpic == Teuchos::null)
      throw(Exceptions::RuntimeError("Cannot cast base Teuchos::Comm to Teuchos::MpiComm object."));
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

    GO ii=0, numAlreadyOwned=0;
    if (myPartitionNumber > -1)
      ii=1;
    //MPI_Reduce(&ii, &numAlreadyOwned, 1, MPI_INT, MPI_SUM, 0,*rawMpiComm);
    sumAll(comm, ii, numAlreadyOwned);
    GetOStream(Statistics0, 0) << "Number of initial arbitration rounds = " << numRounds
                               << ", number partitions assigned = "
                               << numAlreadyOwned << std::endl;

    ///////////////////////////////////////////
    // clean up phase
    ///////////////////////////////////////////

#ifdef MUELU_NEVER_DROP_PROC0
    //If all partitions have been assigned but PID 0 doesn't own one, reassign a partition PID 0.
    if (numAlreadyOwned == numPartitions && mypid==0 && myPartitionNumber==-1 ) {
      // Grab ownership of the partition for which 0 has the most dofs.
      GO maxSize=0;
      for (int i=0; i<allPartitionsIContributeTo.size(); ++i) {
        if (localNnzPerPartition[i] > maxSize) {
          maxSize=localNnzPerPartition[i];
          myPartitionNumber=allPartitionsIContributeTo[i];
        }
      }
      GetOStream(Statistics0, 0) << "Reassigning partition " << myPartitionNumber << " to pid 0." << std::endl;
    }
    // Broadcast PID 0's partition number.  Whichever PID was the previous owner gives up ownership.
    // + It may be that PID 0 didn't own a partition, in which case the number broadcast is -1.
    //   No other PID gives up ownership because partition numbers start at 0.
    // + It may be that PID 0 already owns a partition from the arbitration round.  This is ok
    //   because no other PID will have that partition number.
    int root=0;
    int mpn=-1;
    if (mypid==0) mpn=Teuchos::as<int>(myPartitionNumber);
    MPI_Bcast(&mpn, 1, MPI_INT, root,*rawMpiComm);
    if ( (mypid > 0) && (myPartitionNumber >= 0) && (Teuchos::as<int>(myPartitionNumber) == mpn) ) {
      myPartitionNumber=-1;
    }
#endif //ifdef MUELU_NEVER_DROP_PROC0

    // 1) all pids get global snapshot of partition #s and their owners
    Array<int> allAssignedPartitions(comm->getSize());
    int pn = Teuchos::as<int>(myPartitionNumber);  
    MPI_Allgather(&pn, 1, MPI_INT,
                  allAssignedPartitions.getRawPtr(), 1, MPI_INT,
                  *rawMpiComm);
    //at this point, allAssignedPartitions[i] = partition number owned by pid i

    // 2) Assign unowned partitions to lowest pids who don't already own a partition.
    for (int i=0; i<numPartitions; ++i)
      partitionOwners.push_back(-1);

    int numAdditionalOwners=0;
    std::vector<int> owningPids;
    std::vector<int> eligiblePids;
    eligiblePids.reserve(comm->getSize()-numAlreadyOwned);
    for (int i=0; i<allAssignedPartitions.size(); ++i) {
      if (allAssignedPartitions[i] != -1) {
        partitionOwners[allAssignedPartitions[i]] = i;
      } else { //TODO once we have at least numPartitions-numAlreadyOwned PIDs, we can break
        eligiblePids.push_back(i);
      }
    }

    //if all partitions are already owned, then skip this step.
    if (numAlreadyOwned < numPartitions) {
      int next=0;
      for (int i=0; i<partitionOwners.size(); ++i) {
        if (partitionOwners[i] == -1) {
          partitionOwners[i] = eligiblePids[next++];
          allAssignedPartitions[partitionOwners[i]] = i;
          numAdditionalOwners++;
        }
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION((numAlreadyOwned + numAdditionalOwners) != numPartitions, Exceptions::RuntimeError,
                               "Diffusive clean up phase failed to assign owners to " +
                               Teuchos::toString(numPartitions-numAlreadyOwned-numAdditionalOwners) + " partitions");

    GetOStream(Statistics0, 0) << "Clean up phase: assigned " << numAdditionalOwners << " additional partitions." << std::endl;

    myPartitionNumber = allAssignedPartitions[mypid];

    int numPartitionOwners=0;
    for (int i=0; i<partitionOwners.size(); ++i)
      if (partitionOwners[i] > -1) ++numPartitionOwners;
    TEUCHOS_TEST_FOR_EXCEPTION(numPartitionOwners != numPartitions, Exceptions::RuntimeError,
                               "Number of partition owners (" + Teuchos::toString(numPartitionOwners) + ") is not equal to number of partitions");

    // print the grid of processors
    GetOStream(Statistics1, 0) << "Partition distribution over cores, + indicates partition ownership" << std::endl;
    int numProc = comm->getSize();
    ArrayRCP<char> grid(numProc, '.');
    for (int i=0; i<partitionOwners.size(); ++i) grid[ partitionOwners[i] ] = '+';
    int sizeOfARow = (int) sqrt(numProc);
    int numRows = numProc / sizeOfARow;
    int leftOvers = numProc  - (numProc/sizeOfARow)*sizeOfARow;
    int ctr=0;
    int pidCtr=0;
    for (int i=0; i<numRows; ++i) {
      for (int j=0; j<sizeOfARow; ++j)
        GetOStream(Statistics1, 0) << grid[ctr++];
      GetOStream(Statistics1, 0) << "      " << pidCtr << ":" << pidCtr+sizeOfARow-1 << std::endl;;
      pidCtr += sizeOfARow;
    }
    if (leftOvers > 0) {
      for (int i=0; i<leftOvers; ++i)
        GetOStream(Statistics1, 0) << grid[ctr++];

      Array<char> aos(sizeOfARow-leftOvers, ' ');
      std::string spaces(aos.begin(), aos.end());
      GetOStream(Statistics1, 0) << spaces << "      " << pidCtr << ":" << pidCtr+leftOvers-1 << std::endl;;
    }

  } //DeterminePartitionPlacement

} // namespace MueLu

#endif //ifdef HAVE_MPI

#endif // MUELU_REPARTITIONFACTORY_DEF_HPP
