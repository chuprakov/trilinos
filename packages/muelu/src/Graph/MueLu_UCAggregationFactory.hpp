#ifndef MUELU_UCAGGREGATIONFACTORY_HPP
#define MUELU_UCAGGREGATIONFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_LocalAggregationAlgorithm.hpp"
#include "MueLu_LeftoverAggregationAlgorithm.hpp"

namespace MueLu {

  /*!
    @class UCAggregationFactory class.
    @brief Factory for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.
  */

  /* Factory input:
     - a graph ("Graph")

     Factory output:
     - aggregates ("Aggegates")

     Factory options:
     - 
     - 
     - 
  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class UCAggregationFactory : public SingleLevelFactoryBase {
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    UCAggregationFactory()
    {
      TEST_FOR_EXCEPTION(algo2_.GetMinNodesPerAggregate() != algo1_.GetMinNodesPerAggregate(), Exceptions::RuntimeError, "");
    }

    //! Destructor.
    virtual ~UCAggregationFactory() {}
    //@}

    //! @name Set/get methods.
    //@{

    // Options algo1
    void SetOrdering(Ordering ordering) { algo1_.SetOrdering(ordering); }
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { algo1_.SetMaxNeighAlreadySelected(maxNeighAlreadySelected); }
    Ordering GetOrdering() const { return algo1_.GetOrdering(); }
    int GetMaxNeighAlreadySelected() const { return algo1_.GetmaxNeighAlreadySelected(); }

    // Options algo2
    void SetPhase3AggCreation(double phase3AggCreation) { algo2_.SetPhase3AggCreation(phase3AggCreation); }
    double GetPhase3AggCreation() const { return algo2_.SetPhase3AggCreation(); }

    // Options shared algo1 and algo2
    void SetMinNodesPerAggregate(int minNodesPerAggregate) { algo1_.SetMinNodesPerAggregate(minNodesPerAggregate); algo2_.SetMinNodesPerAggregate(minNodesPerAggregate); }
    int GetMinNodesPerAggregate() const { return algo1_.GetMinNodesPerAggregate(); TEST_FOR_EXCEPTION(algo2_.GetMinNodesPerAggregate() != algo1_.GetMinNodesPerAggregate(), Exceptions::RuntimeError, ""); }
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      //currentLevel.Input("Graph", graphFact_()); // TW: what's this???
      currentLevel.Request("Graph", graphFact_); // we should request data...
    }

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Build aggregates. */
    bool Build(Level &currentLevel) const
    {
      //TODO check for reuse of aggregates here
      //FIXME should there be some way to specify the name of the graph in the needs table, i.e., could
      //FIXME there ever be more than one graph?
      std::ostringstream buf; buf << currentLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("UCAggregationFactory::Build_"+buf.str()));
      timer->start(true);

      // Level Get
      RCP<const Graph> graph = currentLevel.Get< RCP<Graph> >("Graph", graphFact_);

      // Build
      RCP<Aggregates> aggregates = rcp(new Aggregates(*graph, "UC")); 
      algo1_.CoarsenUncoupled(*graph, *aggregates);
      algo2_.AggregateLeftovers(*graph, *aggregates);

      // Level Set
      currentLevel.Set("Aggregates", aggregates, this);

      //
      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(graph->GetComm()));

      return true; // ??
    }

  private:

    //! Graph Factory
    RCP<FactoryBase> graphFact_;
 
    //! Algorithms
    LocalAggregationAlgorithm algo1_;
    LeftoverAggregationAlgorithm algo2_;

  }; // class UCAggregationFactory

} //namespace MueLu

#define MUELU_UCAGGREGATIONFACTORY_SHORT
#endif //ifndef MUELU_UCAGGREGATIONFACTORY_HPP

//TODO: can be more generic:
// - allow to choose algo 
// - base class for algorithm and options forward to algorithm as parameter list
