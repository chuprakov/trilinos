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
 * MueLu_UncoupledAggregationFactory_def.hpp
 *
 *  Created on: Sep 17, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_UncoupledAggregationFactory_decl.hpp"

#include "MueLu_OnePtAggregationAlgorithm.hpp"
#include "MueLu_SmallAggregationAlgorithm.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm.hpp"
#include "MueLu_UncoupledAggregationAlgorithm.hpp"
#include "MueLu_MaxLinkAggregationAlgorithm.hpp"
#include "MueLu_IsolatedNodeAggregationAlgorithm.hpp"
#include "MueLu_EmergencyAggregationAlgorithm.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UncoupledAggregationFactory()
  : bDefinitionPhase_(true)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // input parameters
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");

    // Aggregation parameters (used in aggregation algorithms)
    // TODO introduce local member function for each aggregation algorithm such that each aggregation algorithm can define its own parameters
    validParamList->set<Ordering>("Ordering", AggOptions::NATURAL, "Ordering strategy (NATURAL|GRAPH|RANDOM)");
    validParamList->set<LO>      ("MaxNeighAlreadySelected",    0, "Number of maximum neighbour nodes that are already aggregated already. "
                                  "If a new aggregate has some neighbours that are already aggregated, "
                                  "this node probably can be added to one of these aggregates. We don't need a new one.");
    validParamList->set<LO>      ("MinNodesPerAggregate",       2, "Minimum number of nodes for aggregate");

    validParamList->set<bool> ("UseOnePtAggregationAlgorithm",              true, "Allow special nodes to be marked for one-to-one transfer to the coarsest level. (default = on)");
    validParamList->set<bool> ("UseSmallAggregatesAggregationAlgorithm",   false, "Turn on/off build process for small aggregates in user defined regions. (default = off)");
    validParamList->set<bool> ("UsePreserveDirichletAggregationAlgorithm", false, "Turn on/off aggregate Dirichlet (isolated nodes) into separate 1pt node aggregates (default = off)");
    validParamList->set<bool> ("UseUncoupledAggregationAlgorithm",          true, "Turn on/off uncoupled aggregation process. Do not turn off: this is "
                               "the main aggregation routine within the uncoupled aggregation process. (default = on)");
    validParamList->set<bool> ("UseMaxLinkAggregationAlgorithm",            true, "Turn on/off MaxLink aggregation algorithm. Adds non-aggregated nodes to "
                               "the next already aggregated neighbour node with the most links. (default = on)");
    validParamList->set<bool> ("UseIsolatedNodeAggregationAlgorithm",       true, "Turn on/off IsolatedNode aggregation algorithm. Ignores isolated "
                               "nodes during aggregation process. (default = on)");
    validParamList->set<bool> ("UseEmergencyAggregationAlgorithm",          true, "Turn on/off Emergency aggregation algorithm. Puts all left over nodes "
                               "into aggregates (including very small aggregates or one-point aggregates). (default = on)");

    validParamList->set< std::string >           ("OnePt aggregate map name",                  "", "Name of input map for single node aggregates. (default='')");
    validParamList->set< RCP<const FactoryBase> >("OnePt aggregate map factory",    null, "Generating factory of (DOF) map for single node aggregates.");
    validParamList->set< std::string >           ("SmallAgg aggregate map name",               "", "Name of input map for small aggregates. (default='')");
    validParamList->set< RCP<const FactoryBase> >("SmallAgg aggregate map factory", null, "Generating factory of (DOF) map for small aggregates.");

    return validParamList;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");
    Input(currentLevel, "DofsPerNode");

    const ParameterList& pL = GetParameterList();
    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name"), mapSmallAggName = pL.get<std::string>("SmallAgg aggregate map name");

    if (mapOnePtName.length() > 0) {
      RCP<const FactoryBase> mapOnePtFact = GetFactory("OnePt aggregate map factory");
      currentLevel.DeclareInput(mapOnePtName, mapOnePtFact.get());
    }
    if (mapSmallAggName.length() > 0) {
      RCP<const FactoryBase> mapSmallAggFact = GetFactory("SmallAgg aggregate map factory");
      currentLevel.DeclareInput(mapSmallAggName, mapSmallAggFact.get());
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    const ParameterList& pL = GetParameterList();
    bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

    bool bUseOnePtAggregationAlgorithm             = pL.get<bool>("UseOnePtAggregationAlgorithm");
    bool bUseSmallAggregationAlgorithm             = pL.get<bool>("UseSmallAggregatesAggregationAlgorithm");
    bool bUsePreserveDirichletAggregationAlgorithm = pL.get<bool>("UsePreserveDirichletAggregationAlgorithm");
    bool bUseUncoupledAggregationAglorithm         = pL.get<bool>("UseUncoupledAggregationAlgorithm");
    bool bUseMaxLinkAggregationAlgorithm           = pL.get<bool>("UseMaxLinkAggregationAlgorithm");
    bool bUseIsolatedNodeAggregationAglorithm      = pL.get<bool>("UseIsolatedNodeAggregationAlgorithm");
    bool bUseEmergencyAggregationAlgorithm         = pL.get<bool>("UseEmergencyAggregationAlgorithm");

    // define aggregation algorithms
    RCP<const FactoryBase> graphFact = GetFactory("Graph");

    // TODO Can we keep different aggregation algorithms over more Build calls?
    algos_.clear();
    if (bUseOnePtAggregationAlgorithm)             algos_.push_back(rcp(new OnePtAggregationAlgorithm             (graphFact)));
    if (bUseSmallAggregationAlgorithm)             algos_.push_back(rcp(new SmallAggregationAlgorithm             (graphFact)));
    if (bUseUncoupledAggregationAglorithm)         algos_.push_back(rcp(new UncoupledAggregationAlgorithm         (graphFact)));
    if (bUseMaxLinkAggregationAlgorithm)           algos_.push_back(rcp(new MaxLinkAggregationAlgorithm           (graphFact)));
    if (bUsePreserveDirichletAggregationAlgorithm) algos_.push_back(rcp(new PreserveDirichletAggregationAlgorithm (graphFact)));
    if (bUseIsolatedNodeAggregationAglorithm)      algos_.push_back(rcp(new IsolatedNodeAggregationAlgorithm      (graphFact)));
    if (bUseEmergencyAggregationAlgorithm)         algos_.push_back(rcp(new EmergencyAggregationAlgorithm         (graphFact)));


    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name"), mapSmallAggName = pL.get<std::string>("SmallAgg aggregate map name");
    RCP<const Map> OnePtMap, SmallAggMap;
    if (mapOnePtName.length()) {
      RCP<const FactoryBase> mapOnePtFact = GetFactory("OnePt aggregate map factory");
      OnePtMap = currentLevel.Get<RCP<const Map> >(mapOnePtName, mapOnePtFact.get());
    }
    if (mapSmallAggName.length()) {
      RCP<const FactoryBase> mapSmallAggFact = GetFactory("SmallAgg aggregate map factory");
      SmallAggMap = currentLevel.Get<RCP<const Map> >(mapSmallAggName, mapSmallAggFact.get());
    }

    RCP<const GraphBase> graph = Get< RCP<GraphBase> >(currentLevel, "Graph");

    // Build
    RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("UC");

    const LO nRows = graph->GetNodeNumVertices();

    // construct aggStat information
    std::vector<unsigned> aggStat(nRows, NodeStats::READY);

    ArrayRCP<const bool> dirichletBoundaryMap = graph->GetBoundaryNodeMap();
    if (dirichletBoundaryMap != Teuchos::null) {
      for (LO i = 0; i < nRows; i++)
        if (dirichletBoundaryMap[i] == true)
          aggStat[i] = NodeStats::BOUNDARY;
    }

    LO nDofsPerNode = Get<LO>(currentLevel, "DofsPerNode");
    GO indexBase = graph->GetDomainMap()->getIndexBase();
    if (SmallAggMap != Teuchos::null || OnePtMap != Teuchos::null) {
      for (LO i = 0; i < nRows; i++) {
        // reconstruct global row id (FIXME only works for contiguous maps)
        GO grid = (graph->GetDomainMap()->getGlobalElement(i)-indexBase) * nDofsPerNode + indexBase;

        if (SmallAggMap != null) {
          for (LO kr = 0; kr < nDofsPerNode; kr++) {
            if (SmallAggMap->isNodeGlobalElement(grid + kr))
              aggStat[i] = MueLu::NodeStats::SMALLAGG;
          }
        }

        if (OnePtMap != null) {
          for (LO kr = 0; kr < nDofsPerNode; kr++) {
            if (OnePtMap->isNodeGlobalElement(grid + kr))
              aggStat[i] = MueLu::NodeStats::ONEPT;
          }
        }
      }
    }

    LO numNonAggregatedNodes = nRows;

    for (size_t a = 0; a < algos_.size(); a++) {
      SubFactoryMonitor sfm(*this, "Algo \"" + algos_[a]->description() + "\"", currentLevel);
      algos_[a]->BuildAggregates(pL, *graph, *aggregates, aggStat, numNonAggregatedNodes);

      if (numNonAggregatedNodes == 0) {
        // All nodes have been aggregated, can quit early
        break;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes > 0, Exceptions::RuntimeError, "MueLu::UncoupledAggregationFactory::Build: Leftover nodes found! Error!");

    aggregates->AggregatesCrossProcessors(false);

    Set(currentLevel, "Aggregates", aggregates);

    GetOStream(Statistics0, 0) << aggregates->description() << std::endl;
  }

} //namespace MueLu


#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_ */
