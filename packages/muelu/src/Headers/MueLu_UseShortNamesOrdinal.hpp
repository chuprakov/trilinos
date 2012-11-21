// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Xpetra_UseShortNamesOrdinal.hpp>

#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Graph;
#endif

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Aggregates;
#endif

#ifdef MUELU_LOCALAGGREGATIONALGORITHM_SHORT
typedef MueLu::LocalAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LocalAggregationAlgorithm;
#endif

#ifdef MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT
typedef MueLu::LeftoverAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LeftoverAggregationAlgorithm;
#endif

#ifdef MUELU_UCAGGREGATIONFACTORY_SHORT
typedef MueLu::UCAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UCAggregationFactory;
#endif

#ifdef MUELU_UCAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::UCAggregationCommHelper<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UCAggregationCommHelper;
#endif

#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PRFactory;
#endif

#ifdef MUELU_ZOLTANINTERFACE_SHORT
typedef MueLu::ZoltanInterface<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ZoltanInterface;
#endif

#ifdef MUELU_AMALGAMATIONINFO_SHORT
typedef MueLu::AmalgamationInfo<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AmalgamationInfo;
#endif

#ifdef MUELU_CHEAPAGGREGATIONALGORITHM_SHORT
typedef MueLu::CheapAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CheapAggregationAlgorithm;
#endif

#ifdef MUELU_EXPERIMENTALAGGREGATIONFACTORY_SHORT
typedef MueLu::ExperimentalAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ExperimentalAggregationFactory;
#endif

#ifdef MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
typedef MueLu::UncoupledAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UncoupledAggregationFactory;
#endif

#ifdef MUELU_UNCOUPLEDAGGREGATIONALGORITHM_SHORT
typedef MueLu::UncoupledAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UncoupledAggregationAlgorithm;
#endif

#ifdef MUELU_MAXLINKAGGREGATIONALGORITHM_SHORT
typedef MueLu::MaxLinkAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MaxLinkAggregationAlgorithm;
#endif

#ifdef MUELU_EMERGENCYAGGREGATIONALGORITHM_SHORT
typedef MueLu::EmergencyAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> EmergencyAggregationAlgorithm;
#endif

#ifdef MUELU_ONEPTAGGREGATIONALGORITHM_SHORT
typedef MueLu::OnePtAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> OnePtAggregationAlgorithm;
#endif

#ifdef MUELU_SINGLELEVELFACTORYBASE_SHORT
typedef MueLu::SingleLevelFactoryBase SingleLevelFactoryBase;
#endif

#ifdef MUELU_TWOLEVELFACTORYBASE_SHORT
typedef MueLu::TwoLevelFactoryBase TwoLevelFactoryBase;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory RFactory;
#endif

#ifdef MUELU_AMESOSSMOOTHER_SHORT
typedef MueLu::AmesosSmoother AmesosSmoother;
#endif

#ifdef MUELU_IFPACKSMOOTHER_SHORT
typedef MueLu::IfpackSmoother IfpackSmoother;
#endif

#ifdef MUELU_FACTORYBASE_SHORT
typedef MueLu::FactoryBase FactoryBase;
#endif

#ifdef MUELU_FACTORY_SHORT
typedef MueLu::Factory Factory;
#endif

#ifdef MUELU_FACTORYMANAGERBASE_SHORT
typedef MueLu::FactoryManagerBase FactoryManagerBase;
#endif

#ifdef MUELU_LEVEL_SHORT
typedef MueLu::Level Level;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory PFactory;
#endif

#ifdef MUELU_TWOKEYMAP_SHORT
typedef MueLu::TwoKeyMap TwoKeyMap;
#endif

#ifdef MUELU_VARIABLECONTAINER_SHORT
typedef MueLu::VariableContainer VariableContainer;
#endif

