// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Xpetra_UseShortNamesOrdinal.hpp>

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Aggregates;
#endif

#ifdef MUELU_AGGREGATIONPHASE1ALGORITHM_SHORT
typedef MueLu::AggregationPhase1Algorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AggregationPhase1Algorithm;
#endif

#ifdef MUELU_AGGREGATIONPHASE2AALGORITHM_SHORT
typedef MueLu::AggregationPhase2aAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AggregationPhase2aAlgorithm;
#endif

#ifdef MUELU_AGGREGATIONPHASE2BALGORITHM_SHORT
typedef MueLu::AggregationPhase2bAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AggregationPhase2bAlgorithm;
#endif

#ifdef MUELU_AMALGAMATIONINFO_SHORT
typedef MueLu::AmalgamationInfo<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AmalgamationInfo;
#endif

#ifdef MUELU_COUPLEDAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::CoupledAggregationCommHelper<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoupledAggregationCommHelper;
#endif

#ifdef MUELU_COUPLEDAGGREGATIONFACTORY_SHORT
typedef MueLu::CoupledAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoupledAggregationFactory;
#endif

#ifdef MUELU_DEMOFACTORY_SHORT
typedef MueLu::DemoFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> DemoFactory;
#endif

#ifdef MUELU_EMERGENCYAGGREGATIONALGORITHM_SHORT
typedef MueLu::EmergencyAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> EmergencyAggregationAlgorithm;
#endif

#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Graph;
#endif

#ifdef MUELU_GRAPHBASE_SHORT
typedef MueLu::GraphBase<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> GraphBase;
#endif

#ifdef MUELU_ISOLATEDNODEAGGREGATIONALGORITHM_SHORT
typedef MueLu::IsolatedNodeAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> IsolatedNodeAggregationAlgorithm;
#endif

#ifdef MUELU_ISORROPIAINTERFACE_SHORT
typedef MueLu::IsorropiaInterface<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> IsorropiaInterface;
#endif

#ifdef MUELU_LWGRAPH_SHORT
typedef MueLu::LWGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LWGraph;
#endif

#ifdef MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT
typedef MueLu::LeftoverAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LeftoverAggregationAlgorithm;
#endif

#ifdef MUELU_LOCALAGGREGATIONALGORITHM_SHORT
typedef MueLu::LocalAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LocalAggregationAlgorithm;
#endif

#ifdef MUELU_MAXLINKAGGREGATIONALGORITHM_SHORT
typedef MueLu::MaxLinkAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MaxLinkAggregationAlgorithm;
#endif

#ifdef MUELU_ONEPTAGGREGATIONALGORITHM_SHORT
typedef MueLu::OnePtAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> OnePtAggregationAlgorithm;
#endif

#ifdef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_SHORT
typedef MueLu::PreserveDirichletAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PreserveDirichletAggregationAlgorithm;
#endif

#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PRFactory;
#endif

#ifdef MUELU_REBALANCEMAPFACTORY_SHORT
typedef MueLu::RebalanceMapFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RebalanceMapFactory;
#endif

#ifdef MUELU_REPARTITIONINTERFACE_SHORT
typedef MueLu::RepartitionInterface<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RepartitionInterface;
#endif

#ifdef MUELU_SMALLAGGREGATIONALGORITHM_SHORT
typedef MueLu::SmallAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SmallAggregationAlgorithm;
#endif

#ifdef MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
typedef MueLu::UncoupledAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UncoupledAggregationFactory;
#endif

#ifdef MUELU_USERAGGREGATIONFACTORY_SHORT
typedef MueLu::UserAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UserAggregationFactory;
#endif

#ifdef MUELU_ZOLTAN2INTERFACE_SHORT
typedef MueLu::Zoltan2Interface<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Zoltan2Interface;
#endif

#ifdef MUELU_ZOLTANINTERFACE_SHORT
typedef MueLu::ZoltanInterface<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ZoltanInterface;
#endif

#ifdef MUELU_AMESOSSMOOTHER_SHORT
typedef MueLu::AmesosSmoother AmesosSmoother;
#endif

#ifdef MUELU_FACTORY_SHORT
typedef MueLu::Factory Factory;
#endif

#ifdef MUELU_FACTORYBASE_SHORT
typedef MueLu::FactoryBase FactoryBase;
#endif

#ifdef MUELU_FACTORYMANAGERBASE_SHORT
typedef MueLu::FactoryManagerBase FactoryManagerBase;
#endif

#ifdef MUELU_IFPACKSMOOTHER_SHORT
typedef MueLu::IfpackSmoother IfpackSmoother;
#endif

#ifdef MUELU_LEVEL_SHORT
typedef MueLu::Level Level;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory PFactory;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory RFactory;
#endif

#ifdef MUELU_SINGLELEVELFACTORYBASE_SHORT
typedef MueLu::SingleLevelFactoryBase SingleLevelFactoryBase;
#endif

#ifdef MUELU_TWOLEVELFACTORYBASE_SHORT
typedef MueLu::TwoLevelFactoryBase TwoLevelFactoryBase;
#endif

#ifdef MUELU_VARIABLECONTAINER_SHORT
typedef MueLu::VariableContainer VariableContainer;
#endif

#ifdef MUELU_SMOOTHERFACTORYBASE_SHORT
typedef MueLu::SmootherFactoryBase SmootherFactoryBase;
#endif
