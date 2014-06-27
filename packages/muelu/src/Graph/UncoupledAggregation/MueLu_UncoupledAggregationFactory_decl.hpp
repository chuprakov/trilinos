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
#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_DECL_HPP_


#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"

#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_OnePtAggregationAlgorithm_fwd.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm_fwd.hpp"
#include "MueLu_MaxLinkAggregationAlgorithm_fwd.hpp"
#include "MueLu_IsolatedNodeAggregationAlgorithm_fwd.hpp"
#include "MueLu_EmergencyAggregationAlgorithm_fwd.hpp"

#include "MueLu_AggregationPhase1Algorithm_fwd.hpp"
#include "MueLu_AggregationPhase2aAlgorithm_fwd.hpp"
#include "MueLu_AggregationPhase2bAlgorithm_fwd.hpp"
#include "MueLu_AggregationPhase3Algorithm_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
//#include "MueLu_Graph_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType, class LocalMatOps = typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class UncoupledAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  UncoupledAggregationFactory();

  //! Destructor.
  virtual ~UncoupledAggregationFactory() { }

  RCP<const ParameterList> GetValidParameterList(const ParameterList& paramList = ParameterList()) const;

  //@}

  //! @name Set/get methods.
  //@{

  // Options shared by all aggregation algorithms

  // deprecated
  void SetOrdering(AggOptions::Ordering ordering) {
    SetParameter("Ordering", ParameterEntry(ordering)); // revalidate
  }
  // deprecated
  void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) {
    SetParameter("MaxNeighAlreadySelected", ParameterEntry(Teuchos::as<LocalOrdinal>(maxNeighAlreadySelected))); // revalidate
  }
  // deprecated
  void SetMinNodesPerAggregate(int minNodesPerAggregate) {
    SetParameter("MinNodesPerAggregate", ParameterEntry(Teuchos::as<LocalOrdinal>(minNodesPerAggregate))); // revalidate
  }
  // set information about 1-node aggregates (map name and generating factory)
  void SetOnePtMapName(const std::string name, Teuchos::RCP<const FactoryBase> mapFact) {
    SetParameter("OnePt aggregate map name", ParameterEntry(std::string(name))); // revalidate
    SetFactory("OnePt aggregate map factory",mapFact);
  }

  // deprecated
  Ordering GetOrdering() const {
    const ParameterList& pL = GetParameterList();
    Ordering ordering = pL.get<Ordering>("Ordering");
    return ordering;
  }
  // deprecated
  int GetMaxNeighAlreadySelected() const {
    const ParameterList& pL = GetParameterList();
    return Teuchos::as<int>(pL.get<LocalOrdinal>("MaxNeighAlreadySelected"));
  }
  // deprecated
  int GetMinNodesPerAggregate() const {
    const ParameterList& pL = GetParameterList();
    return Teuchos::as<int>(pL.get<LocalOrdinal>("MinNodesPerAggregate"));
  }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level &currentLevel) const;

  //@}

  //! @name Definition methods
  //@{

  /*! @brief Append a new aggregation algorithm to list of aggregation algorithms */
  //void Append(const RCP<MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & alg);

  /*! @brief Remove all aggregation algorithms from list */
  //void ClearAggregationAlgorithms() { algos_.clear(); }
  //@}

private:

  //! aggregation algorithms
  // will be filled in Build routine
  mutable std::vector<RCP<MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > > algos_;

  //! boolean flag: definition phase
  //! if true, the aggregation algorithms still can be set and changed.
  //! if false, no change in aggregation algorithms is possible any more
  mutable bool bDefinitionPhase_;

}; // class UncoupledAggregationFactory

}

#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DECL_HPP_ */
