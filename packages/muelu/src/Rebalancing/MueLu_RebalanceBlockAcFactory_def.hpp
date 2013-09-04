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
 * MueLu_RebalanceBlockAcFactory_def.hpp
 *
 *  Created on: Aug 15, 2013
 *      Author: tobias
 */

#ifndef MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_RAPFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
#include "MueLu_RepartitionFactory.hpp"
//#include "MueLu_RebalanceTransferFactory.hpp"
//#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_IsorropiaInterface.hpp"
//#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#endif

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RebalanceBlockAcFactory() {  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",         Teuchos::null, "Generating factory of the matrix A for rebalancing");
    //validParamList->set< bool >                  ("useSubcomm",         true, "Construct subcommunicators");
    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
    FactManager_.push_back(FactManager);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(coarseLevel, "A");

    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
      SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

      coarseLevel.DeclareInput("Importer",(*it)->GetFactory("Importer").get(), this);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Computing blocked Ac", coarseLevel);

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    RCP<Matrix> originalAc = Get< RCP<Matrix> >(coarseLevel, "A");

    RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > bA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(originalAc);
    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

    // plausibility check
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != 2,Exceptions::RuntimeError, "MueLu::BlockedPFactory::Build: number of block rows of A is not equal 2. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols() != 2,Exceptions::RuntimeError, "MueLu::BlockedPFactory::Build: number of block columns of A is not equal 2. error.");

    // store map extractors
    Teuchos::RCP<const MapExtractorClass> rangeMapExtractor  = bA->getRangeMapExtractor();
    Teuchos::RCP<const MapExtractorClass> domainMapExtractor = bA->getDomainMapExtractor();

    std::vector<GO> fullRangeMapVector;
    std::vector<GO> fullDomainMapVector;

    std::vector<RCP<const Map> > subBlockARangeMaps;
    std::vector<RCP<const Map> > subBlockADomainMaps;
    subBlockARangeMaps.reserve(bA->Rows());
    subBlockADomainMaps.reserve(bA->Cols());

    std::vector<Teuchos::RCP<Matrix> > subBlockRebA;
    subBlockRebA.reserve(bA->Cols() * bA->Rows());

    for(size_t i=0; i<bA->Rows(); i++) {
      for(size_t j=0; j<bA->Cols(); j++) {
        // extract matrix block
        Teuchos::RCP<CrsMatrix> Amij = bA->getMatrix(i, j);
        Teuchos::RCP<CrsMatrixWrap> Awij = Teuchos::rcp(new CrsMatrixWrap(Amij));
        Teuchos::RCP<Matrix> Aij = Teuchos::rcp_dynamic_cast<Matrix>(Awij);
        //subBlockRebA[i*bA->Cols() + j] = Aij;
        subBlockRebA.push_back(Aij);
      }
    }

    size_t curBlockId = 0;
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
      SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

      Teuchos::RCP<const Import> rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());

      // rebalance diagonal block

      // extract matrix block
      Teuchos::RCP<Matrix> Aii = subBlockRebA[curBlockId*bA->Cols() + curBlockId];

      Teuchos::RCP<Matrix> rebAii;
      if(rebalanceImporter != Teuchos::null) {
        std::stringstream ss; ss << "Rebalancing matrix block A(" << curBlockId << "," << curBlockId << ")";
        SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
        RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

        //const ParameterList & pL = GetParameterList();

        ParameterList XpetraList;
        //if (pL.get<bool>("useSubcomm") == true) {
          //GetOStream(Runtime0,0) << "Replacing maps with a subcommunicator" << std::endl;
          XpetraList.set("Restrict Communicator",false /*true*/ /*XXX*/);
        //}
        // NOTE: If the communicator is restricted away, Build returns Teuchos::null.
        rebAii = MatrixFactory::Build(Aii, *rebalanceImporter, targetMap, targetMap, rcp(&XpetraList,false));

        if (!rebAii.is_null())
          rebAii->SetFixedBlockSize(Aii->GetFixedBlockSize());

        if (!rebAii.is_null()) {
          RCP<ParameterList> params = rcp(new ParameterList());
          params->set("printLoadBalancingInfo", true);
          std::stringstream ss2; ss2 << "A(" << curBlockId << "," << curBlockId << ") rebalanced:";
          GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAii, ss2.str(), params);
        }
      }  // rebalance matrix block A(i,i)
      else {
        rebAii = Aii;
        /*RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        std::stringstream ss2; ss2 << "A(" << curBlockId << "," << curBlockId << ") not rebalanced:";
        GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAii, ss2.str(), params);*/
      }

      subBlockRebA[curBlockId*bA->Cols() + curBlockId] = rebAii;

      // rebalance off-diagonal matrix blocks in same row
      for(size_t j=0; j<bA->Cols(); j++) {
        if(j==curBlockId) continue;  // jump over block diagonal matrix block

        // extract matrix block
        Teuchos::RCP<Matrix> Aij = subBlockRebA[curBlockId*bA->Cols() + j];

        Teuchos::RCP<Matrix> rebAij;
        if(rebalanceImporter!=Teuchos::null) {
          std::stringstream ss3; ss3 << "Rebalancing matrix block A(" << curBlockId << "," << j << ")";
          SubFactoryMonitor subM(*this, ss3.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          Teuchos::RCP<Map> dummy; // The 3rd argument says to use the original domain map
          rebAij = MatrixFactory::Build(Aij, *rebalanceImporter, dummy, targetMap);

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebAij.is_null() && Aij->IsView("stridedMaps"))
            rebAij->CreateView("stridedMaps", Aij);

          if (!rebAij.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            std::stringstream ss4; ss4 << "A(" << curBlockId << "," << j << ") rebalanced:";
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss4.str(), params);
          }
        } // rebalance matrix block A(i,j)
        else {
          rebAij = Aij;
          /*RCP<ParameterList> params = rcp(new ParameterList());
          params->set("printLoadBalancingInfo", true);
          std::stringstream ss2; ss2 << "A(" << curBlockId << "," << j << ") not rebalanced:";
          GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss2.str(), params);*/
        }

        subBlockRebA[curBlockId*bA->Cols() + j] = rebAij;
      } // end loop over all columns

      // rebalance off-diagonal matrix blocks in same column
      for(size_t i=0; i<bA->Rows(); i++) {
        if(i==curBlockId) continue;  // jump over block diagonal matrix block

        // extract matrix block
        Teuchos::RCP<Matrix> Aij = subBlockRebA[i*bA->Cols() + curBlockId];

        Teuchos::RCP<Matrix> rebAij;
        if(rebalanceImporter!=Teuchos::null) {
          std::stringstream ss; ss << "Rebalancing matrix block (" << i << "," << curBlockId << ")";
          SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          rebAij = Aij; // just a copy
          Teuchos::RCP<const CrsMatrixWrap> rebAwij = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebAij);
          Teuchos::RCP<CrsMatrix> rebAmij = rebAwij->getCrsMatrix();
          Teuchos::RCP<const Import> rebAijImport = ImportFactory::Build(targetMap,Aij->getColMap());
          rebAmij->replaceDomainMapAndImporter(targetMap,rebAijImport);

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebAij.is_null() && Aij->IsView("stridedMaps"))
            rebAij->CreateView("stridedMaps", Aij);

          if (!rebAij.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            std::stringstream ss2; ss2 << "A(" << i << "," << curBlockId << ") rebalanced:";
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss2.str(), params);
          }
        } // rebalance matrix block A(1,0)
        else {
          rebAij = Aij;
          /*RCP<ParameterList> params = rcp(new ParameterList());
          params->set("printLoadBalancingInfo", true);
          std::stringstream ss2; ss2 << "A(" << i << "," << curBlockId << ") not rebalanced:";
          GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss2.str(), params);*/
        }

        subBlockRebA[i*bA->Cols() + curBlockId] = rebAij;
      } // end loop over all rows


      // build full range and full domain map (strided)
      subBlockARangeMaps.push_back(rebAii->getRangeMap());
      Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap = subBlockARangeMaps.back()->getNodeElementList();
      fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap.begin(), nodeRangeMap.end());
      sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

      subBlockADomainMaps.push_back(rebAii->getDomainMap());
      Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap = subBlockADomainMaps.back()->getNodeElementList();
      fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap.begin(), nodeDomainMap.end());
      sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

      curBlockId++;
    } // end loop over all block rows

    // now, subBlockRebA contains all rebalanced matrix blocks

    // extract map index base from maps of blocked A
    GO rangeIndexBase  = 0;
    GO domainIndexBase = 0;
    rangeIndexBase = bA->getRangeMap()->getIndexBase();
    domainIndexBase= bA->getDomainMap()->getIndexBase();

    Teuchos::ArrayView<GO> fullRangeMapGIDs(&fullRangeMapVector[0],fullRangeMapVector.size());
    Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getFullMap());
    std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
    Teuchos::RCP<const StridedMap > fullRangeMap =
        StridedMapFactory::Build(
            bA->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            stridedData,
            bA->getRangeMap()->getComm(),
            stridedRgFullMap->getStridedBlockId(),
            stridedRgFullMap->getOffset());

    Teuchos::ArrayView<GO> fullDomainMapGIDs(&fullDomainMapVector[0],fullDomainMapVector.size());

    Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getFullMap());
    std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
    Teuchos::RCP<const StridedMap > fullDomainMap =
        StridedMapFactory::Build(
            bA->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            stridedData2,
            bA->getDomainMap()->getComm(),
            stridedDoFullMap->getStridedBlockId(),
            stridedDoFullMap->getOffset());

    // build map extractors
    Teuchos::RCP<const MapExtractorClass> rebRangeMapExtractor  = MapExtractorFactoryClass::Build(fullRangeMap, subBlockARangeMaps);
    Teuchos::RCP<const MapExtractorClass> rebDomainMapExtractor = MapExtractorFactoryClass::Build(fullDomainMap, subBlockADomainMaps);

    Teuchos::RCP<BlockedCrsMatrix> reb_bA = Teuchos::rcp(new BlockedCrsMatrix(rebRangeMapExtractor,rebDomainMapExtractor,10));
    for(size_t i=0; i<bA->Rows(); i++) {
      for(size_t j=0; j<bA->Cols(); j++) {
       Teuchos::RCP<const CrsMatrixWrap> crsOpij = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(subBlockRebA[i*bA->Cols() + j]);
       Teuchos::RCP<CrsMatrix> crsMatij = crsOpij->getCrsMatrix();
       reb_bA->setMatrix(i,j,crsMatij);
      }
    }
    reb_bA->fillComplete();
    //reb_bA->describe(*out,Teuchos::VERB_EXTREME);
    coarseLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(reb_bA), this);


#if 0
    if (rebalanceFacts_.begin() != rebalanceFacts_.end()) {
      SubFactoryMonitor m2(*this, "Rebalance additional data", coarseLevel);

      // call Build of all user-given transfer factories
      for (std::vector<RCP<const FactoryBase> >::const_iterator it = rebalanceFacts_.begin(); it != rebalanceFacts_.end(); ++it) {
        GetOStream(Runtime0, 0) << "RebalanceAc: call rebalance factory " << (*it).get() << ": " << (*it)->description() << std::endl;
        (*it)->CallBuild(coarseLevel);
      }
    }
#endif
  } //Build()

#if 0
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockedAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddRebalanceFactory(const RCP<const FactoryBase>& factory) {

    /*TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                               "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                               "This is very strange. (Note: you can remove this exception if there's a good reason for)");
    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");*/
    rebalanceFacts_.push_back(factory);
  } //AddRebalanceFactory()
#endif
} //namespace MueLu


#endif /* MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_ */
