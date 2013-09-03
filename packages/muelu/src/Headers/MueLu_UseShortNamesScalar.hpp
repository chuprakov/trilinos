// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

#include <Xpetra_UseShortNamesScalar.hpp>

#ifdef MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_SHORT
typedef MueLu::AdaptiveSaMLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AdaptiveSaMLParameterListInterpreter;
#endif

#ifdef MUELU_AGGREGATIONEXPORTFACTORY_SHORT
typedef MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AggregationExportFactory;
#endif

#ifdef MUELU_AMALGAMATIONFACTORY_SHORT
typedef MueLu::AmalgamationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AmalgamationFactory;
#endif

#ifdef MUELU_AMESOS2SMOOTHER_SHORT
typedef MueLu::Amesos2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Amesos2Smoother;
#endif

#ifdef MUELU_ALGEBRAICPERMUTATIONSTRATEGY_SHORT
typedef MueLu::AlgebraicPermutationStrategy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AlgebraicPermutationStrategy;
#endif

#ifdef MUELU_BLOCKEDCOARSEMAPFACTORY_SHORT
typedef MueLu::BlockedCoarseMapFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> BlockedCoarseMapFactory;
#endif

#ifdef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_SHORT
typedef MueLu::BlockedGaussSeidelSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> BlockedGaussSeidelSmoother;
#endif

#ifdef MUELU_BLOCKEDPFACTORY_SHORT
typedef MueLu::BlockedPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> BlockedPFactory;
#endif

#ifdef MUELU_BLOCKEDRAPFACTORY_SHORT
typedef MueLu::BlockedRAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> BlockedRAPFactory;
#endif

#ifdef MUELU_BRICKAGGREGATIONFACTORY_SHORT
typedef MueLu::BrickAggregationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> BrickAggregationFactory;
#endif

#ifdef MUELU_BRAESSSARAZINSMOOTHER_SHORT
typedef MueLu::BraessSarazinSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> BraessSarazinSmoother;
#endif

#ifdef MUELU_CGSOLVER_SHORT
typedef MueLu::CGSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CGSolver;
#endif

#ifdef MUELU_COALESCEDROPFACTORY_SHORT
typedef MueLu::CoalesceDropFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoalesceDropFactory;
#endif

#ifdef MUELU_COARSEMAPFACTORY_SHORT
typedef MueLu::CoarseMapFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoarseMapFactory;
#endif

#ifdef MUELU_CONSTRAINT_SHORT
typedef MueLu::Constraint<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Constraint;
#endif

#ifdef MUELU_CONSTRAINTFACTORY_SHORT
typedef MueLu::ConstraintFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ConstraintFactory;
#endif

#ifdef MUELU_COORDINATESTRANSFERFACTORY_SHORT
typedef MueLu::CoordinatesTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoordinatesTransferFactory;
#endif

#ifdef MUELU_COUPLEDRBMFACTORY_SHORT
typedef MueLu::CoupledRBMFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoupledRBMFactory;
#endif

#ifdef MUELU_DEMOFACTORY_SHORT
typedef MueLu::DemoFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> DemoFactory;
#endif

#ifdef MUELU_DIRECTSOLVER_SHORT
typedef MueLu::DirectSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> DirectSolver;
#endif

#ifdef MUELU_DISTANCELAPLACIANFACTORY_SHORT
typedef MueLu::DistanceLaplacianFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> DistanceLaplacianFactory;
#endif

#ifdef MUELU_EMINPFACTORY_SHORT
typedef MueLu::EminPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> EminPFactory;
#endif

#ifdef MUELU_FACTORYFACTORY_SHORT
typedef MueLu::FactoryFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> FactoryFactory;
#endif

#ifdef MUELU_FACTORYMANAGER_SHORT
typedef MueLu::FactoryManager<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> FactoryManager;
#endif

#ifdef MUELU_FAKESMOOTHERPROTOTYPE_SHORT
typedef MueLu::FakeSmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> FakeSmootherPrototype;
#endif

#ifdef MUELU_FILTEREDAFACTORY_SHORT
typedef MueLu::FilteredAFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> FilteredAFactory;
#endif

#ifdef MUELU_GENERICRFACTORY_SHORT
typedef MueLu::GenericRFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> GenericRFactory;
#endif

#ifdef MUELU_HIERARCHY_SHORT
typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Hierarchy;
#endif

#ifdef MUELU_HIERARCHYMANAGER_SHORT
typedef MueLu::HierarchyManager<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> HierarchyManager;
#endif

#ifdef MUELU_HIERARCHYFACTORY_SHORT
typedef MueLu::HierarchyFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> HierarchyFactory;
#endif

#ifdef MUELU_IFPACK2SMOOTHER_SHORT
typedef MueLu::Ifpack2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Ifpack2Smoother;
#endif

#ifdef MUELU_LOCALPERMUTATIONSTRATEGY_SHORT
typedef MueLu::LocalPermutationStrategy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LocalPermutationStrategy;
#endif

#ifdef MUELU_MAPTRANSFERFACTORY_SHORT
typedef MueLu::MapTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MapTransferFactory;
#endif

#ifdef MUELU_MERGEDSMOOTHER_SHORT
typedef MueLu::MergedSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MergedSmoother;
#endif

#ifdef MUELU_MLPARAMETERLISTINTERPRETER_SHORT
typedef MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MLParameterListInterpreter;
#endif

#ifdef MUELU_MULTIVECTORTRANSFERFACTORY_SHORT
typedef MueLu::MultiVectorTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MultiVectorTransferFactory;
#endif

#ifdef MUELU_NULLSPACEFACTORY_SHORT
typedef MueLu::NullspaceFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> NullspaceFactory;
#endif

#ifdef MUELU_NULLSPACEPRESMOOTHFACTORY_SHORT
typedef MueLu::NullspacePresmoothFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> NullspacePresmoothFactory;
#endif

#ifdef MUELU_PARAMETERLISTINTERPRETER_SHORT
typedef MueLu::ParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ParameterListInterpreter;
#endif

#ifdef MUELU_PATTERNFACTORY_SHORT
typedef MueLu::PatternFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PatternFactory;
#endif

#ifdef MUELU_PERMUTATIONFACTORY_SHORT
typedef MueLu::PermutationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PermutationFactory;
#endif

#ifdef MUELU_PERMUTINGSMOOTHER_SHORT
typedef MueLu::PermutingSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PermutingSmoother;
#endif

#ifdef MUELU_PGPFACTORY_SHORT
typedef MueLu::PgPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PgPFactory;
#endif

#ifdef MUELU_PREDROPFUNCTIONBASECLASS_SHORT
typedef MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PreDropFunctionBaseClass;
#endif

#ifdef MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
typedef MueLu::PreDropFunctionConstVal<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PreDropFunctionConstVal;
#endif

#ifdef MUELU_PROJECTORSMOOTHER_SHORT
typedef MueLu::ProjectorSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ProjectorSmoother;
#endif

#ifdef MUELU_RAPFACTORY_SHORT
typedef MueLu::RAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RAPFactory;
#endif

#ifdef MUELU_RAPSHIFTFACTORY_SHORT
typedef MueLu::RAPShiftFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RAPShiftFactory;
#endif

#ifdef MUELU_REBALANCEACFACTORY_SHORT
typedef MueLu::RebalanceAcFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RebalanceAcFactory;
#endif

#ifdef MUELU_REBALANCEBLOCKACFACTORY_SHORT
typedef MueLu::RebalanceBlockAcFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RebalanceBlockAcFactory;
#endif

#ifdef MUELU_REBALANCEBLOCKTRANSFERFACTORY_SHORT
typedef MueLu::RebalanceBlockTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RebalanceBlockTransferFactory;
#endif

#ifdef MUELU_REBALANCETRANSFERFACTORY_SHORT
typedef MueLu::RebalanceTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RebalanceTransferFactory;
#endif

#ifdef MUELU_REPARTITIONFACTORY_SHORT
typedef MueLu::RepartitionFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RepartitionFactory;
#endif

#ifdef MUELU_RIGIDBODYMODEFACTORY_SHORT
typedef MueLu::RigidBodyModeFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> RigidBodyModeFactory;
#endif

#ifdef MUELU_SAPFACTORY_SHORT
typedef MueLu::SaPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SaPFactory;
#endif

#ifdef MUELU_SCHURCOMPLEMENTFACTORY_SHORT
typedef MueLu::SchurComplementFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SchurComplementFactory;
#endif

#ifdef MUELU_SHIFTEDLAPLACIAN_SHORT
typedef MueLu::ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ShiftedLaplacian;
#endif

#ifdef MUELU_SHIFTEDLAPLACIANOPERATOR_SHORT
typedef MueLu::ShiftedLaplacianOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ShiftedLaplacianOperator;
#endif

#ifdef MUELU_SIMPLESMOOTHER_SHORT
typedef MueLu::SimpleSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SimpleSmoother;
#endif

#ifdef MUELU_SMOOTHER_SHORT
typedef MueLu::Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Smoother;
#endif

#ifdef MUELU_SMOOTHERBASE_SHORT
typedef MueLu::SmootherBase<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SmootherBase;
#endif

#ifdef MUELU_SMOOTHERFACTORY_SHORT
typedef MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SmootherFactory;
#endif

#ifdef MUELU_SMOOTHERPROTOTYPE_SHORT
typedef MueLu::SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SmootherPrototype;
#endif

#ifdef MUELU_SOLVERBASE_SHORT
typedef MueLu::SolverBase<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SolverBase;
#endif

#ifdef MUELU_STEEPESTDESCENTSOLVER_SHORT
typedef MueLu::SteepestDescentSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SteepestDescentSolver;
#endif

#ifdef MUELU_SUBBLOCKAFACTORY_SHORT
typedef MueLu::SubBlockAFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> SubBlockAFactory;
#endif

#ifdef MUELU_TENTATIVEPFACTORY_SHORT
typedef MueLu::TentativePFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TentativePFactory;
#endif

#ifdef MUELU_THRESHOLDAFILTERFACTORY_SHORT
typedef MueLu::ThresholdAFilterFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ThresholdAFilterFactory;
#endif

#ifdef MUELU_TPETRAOPERATOR_SHORT
typedef MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TpetraOperator;
#endif

#ifdef MUELU_TRANSPFACTORY_SHORT
typedef MueLu::TransPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TransPFactory;
#endif

#ifdef MUELU_TRILINOSSMOOTHER_SHORT
typedef MueLu::TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TrilinosSmoother;
#endif

#ifdef MUELU_USERPFACTORY_SHORT
typedef MueLu::UserPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UserPFactory;
#endif

#ifdef MUELU_UTILITIES_SHORT
typedef MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Utils;
typedef MueLu::Utils2<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Utils2;
#endif
