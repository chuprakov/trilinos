multigrid algorithm = emin
emin: iterative method = sd
emin: num iterations = 4
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
emin: pattern = AkPtent   [default]
problem: symmetric = 1   [default]
aggregation: visualize = 0   [default]
repartition: enable = 0   [default]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 1   [unused]
 relaxation: damping factor = 1   [unused]
 
Level 1
 Prolongator minimization (MueLu::EminPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]
    
   aggregation: mode = old   [default]
   aggregation: max agg size = 2147483647   [default]
   aggregation: min agg size = 2   [default]
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   UseOnePtAggregationAlgorithm = 0   [default]
   UsePreserveDirichletAggregationAlgorithm = 0   [default]
   UseUncoupledAggregationAlgorithm = 1   [default]
   UseMaxLinkAggregationAlgorithm = 1   [default]
   UseIsolatedNodeAggregationAlgorithm = 1   [default]
   UseEmergencyAggregationAlgorithm = 1   [default]
   aggregation: preserve Dirichlet points = 0   [default]
   aggregation: enable phase 1 = 1   [default]
   aggregation: enable phase 2a = 1   [default]
   aggregation: enable phase 2b = 1   [default]
   aggregation: enable phase 3 = 1   [default]
   OnePt aggregate map name =    [default]
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
  Constraint (MueLu::ConstraintFactory)
   Ppattern (MueLu::PatternFactory)
   emin: pattern order = 1   [default]
   
  [empty list]
  
 emin: num iterations = 4
 emin: iterative method = sd
 emin: num reuse iterations = 1   [default]
 P0 = Teuchos::RCP<Xpetra::Matrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [default]
 Keep P0 = 0   [default]
 Constraint0 = Teuchos::RCP<MueLu::Constraint<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [default]
 Keep Constraint0 = 0   [default]
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 1   [unused]
 relaxation: damping factor = 1   [unused]
 
Level 2
 Prolongator minimization (MueLu::EminPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]
    
   aggregation: mode = old   [default]
   aggregation: max agg size = 2147483647   [default]
   aggregation: min agg size = 2   [default]
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   UseOnePtAggregationAlgorithm = 0   [default]
   UsePreserveDirichletAggregationAlgorithm = 0   [default]
   UseUncoupledAggregationAlgorithm = 1   [default]
   UseMaxLinkAggregationAlgorithm = 1   [default]
   UseIsolatedNodeAggregationAlgorithm = 1   [default]
   UseEmergencyAggregationAlgorithm = 1   [default]
   aggregation: preserve Dirichlet points = 0   [default]
   aggregation: enable phase 1 = 1   [default]
   aggregation: enable phase 2a = 1   [default]
   aggregation: enable phase 2b = 1   [default]
   aggregation: enable phase 3 = 1   [default]
   OnePt aggregate map name =    [default]
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
  Constraint (MueLu::ConstraintFactory)
   Ppattern (MueLu::PatternFactory)
   emin: pattern order = 1   [default]
   
  [empty list]
  
 emin: num iterations = 4
 emin: iterative method = sd
 emin: num reuse iterations = 1   [default]
 P0 = Teuchos::RCP<Xpetra::Matrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [default]
 Keep P0 = 0   [default]
 Constraint0 = Teuchos::RCP<MueLu::Constraint<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [default]
 Keep Constraint0 = 0   [default]
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
 Setup Smoother (MueLu::AmesosSmoother{type = Superlu})
 presmoother -> 
  [empty list]
 
 
 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 3
 Operator complexity = 1.44
 Max Coarse Size     = 2000
 Implicit Transpose  = false
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1
 A 1    3333   9997     3.00  1
 A 2    1111   3331     3.00  1
 
 Smoother (level 0) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 1) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 2) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 2) post : no smoother
 
