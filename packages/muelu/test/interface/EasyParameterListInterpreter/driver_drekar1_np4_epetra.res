verbosity = test
number of equations = 1
coarse: max size = 1000
sa: use filtered matrix = 1
filtered matrix: use lumping = 1
smoother: type = CHEBYSHEV
aggregation: drop scheme = laplacian
repartition: enable = 1
repartition: min rows per proc = 2000
repartition: max imbalance = 1.327
repartition: start level = 1
repartition: remap parts = 1
repartition: keep proc 0 = 1
repartition: partitioner = zoltan2
max levels = 10   [default]
debug: graph level = -1   [default]
repartition: rebalance P and R = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
problem: symmetric = 1   [default]
aggregation: visualize = 0   [default]
smoother: params -> 
 chebyshev: degree = 2   [unused]
 chebyshev: ratio eigenvalue = 20   [unused]
 chebyshev: min eigenvalue = 1   [unused]
 chebyshev: zero starting solution = 1   [unused]
 chebyshev: eigenvalue max iterations = 10   [unused]
repartition: params -> 
 algorithm = multijagged   [unused]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 
Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      aggregation threshold = 0
      Dirichlet detection threshold = 0
      algorithm = laplacian
      
     lumping = 1   [unused]
     filtered matrix: reuse graph = 1   [unused]
     filtered matrix: reuse eigenvalue = 1   [unused]
     
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      mode = old   [unused]
      Ordering = 0   [unused]
      MaxNeighAlreadySelected = 0   [unused]
      MinNodesPerAggregate = 2   [unused]
      MaxNodesPerAggregate = 2147483647   [unused]
      UseOnePtAggregationAlgorithm = 0   [unused]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [unused]
      UseMaxLinkAggregationAlgorithm = 1   [unused]
      UseIsolatedNodeAggregationAlgorithm = 1   [unused]
      UseEmergencyAggregationAlgorithm = 1   [unused]
      aggregation: preserve Dirichlet points = 0   [unused]
      aggregation: enable phase 1 = 1   [unused]
      aggregation: enable phase 2a = 1   [unused]
      aggregation: enable phase 2b = 1   [unused]
      aggregation: enable phase 3 = 1   [unused]
      OnePt aggregate map name = 
      
      Build (MueLu::AmalgamationFactory)
      [empty list]
      
      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace
      
      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]
      
     [empty list]
     
    Damping factor = 1.33333
    
    Transpose P (MueLu::TransPFactory)
    [empty list]
    
   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]
   
   Keep AP Pattern = 0
   Keep RAP Pattern = 0
   implicit transpose = 0
   CheckMainDiagonal = 0
   RepairMainDiagonal = 0
   
  startLevel = 1
  minRowsPerProcessor = 2000
  nonzeroImbalance = 1.327
  remapPartitions = 1
  numRemapValues = 4
  alwaysKeepProc0 = 1
  repartition: print partition distribution = 0
  
 type = Interpolation
 implicit = 0
 implicit transpose = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 implicit = 0
 implicit transpose = 0
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20   [unused]
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 
Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      aggregation threshold = 0
      Dirichlet detection threshold = 0
      algorithm = laplacian
      
     lumping = 1   [unused]
     filtered matrix: reuse graph = 1   [unused]
     filtered matrix: reuse eigenvalue = 1   [unused]
     
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      mode = old   [unused]
      Ordering = 0   [unused]
      MaxNeighAlreadySelected = 0   [unused]
      MinNodesPerAggregate = 2   [unused]
      MaxNodesPerAggregate = 2147483647   [unused]
      UseOnePtAggregationAlgorithm = 0   [unused]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [unused]
      UseMaxLinkAggregationAlgorithm = 1   [unused]
      UseIsolatedNodeAggregationAlgorithm = 1   [unused]
      UseEmergencyAggregationAlgorithm = 1   [unused]
      aggregation: preserve Dirichlet points = 0   [unused]
      aggregation: enable phase 1 = 1   [unused]
      aggregation: enable phase 2a = 1   [unused]
      aggregation: enable phase 2b = 1   [unused]
      aggregation: enable phase 3 = 1   [unused]
      OnePt aggregate map name = 
      
      Build (MueLu::AmalgamationFactory)
      [empty list]
      
      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace
      
      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]
      
     [empty list]
     
    Damping factor = 1.33333
    
    Transpose P (MueLu::TransPFactory)
    [empty list]
    
   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]
   
   Keep AP Pattern = 0
   Keep RAP Pattern = 0
   implicit transpose = 0
   CheckMainDiagonal = 0
   RepairMainDiagonal = 0
   
  startLevel = 1
  minRowsPerProcessor = 2000
  nonzeroImbalance = 1.327
  remapPartitions = 1
  numRemapValues = 4
  alwaysKeepProc0 = 1
  repartition: print partition distribution = 0
  
 type = Interpolation
 implicit = 0
 implicit transpose = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 implicit = 0
 implicit transpose = 0
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20   [unused]
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 
Level 3
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      aggregation threshold = 0
      Dirichlet detection threshold = 0
      algorithm = laplacian
      
     lumping = 1   [unused]
     filtered matrix: reuse graph = 1   [unused]
     filtered matrix: reuse eigenvalue = 1   [unused]
     
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      mode = old   [unused]
      Ordering = 0   [unused]
      MaxNeighAlreadySelected = 0   [unused]
      MinNodesPerAggregate = 2   [unused]
      MaxNodesPerAggregate = 2147483647   [unused]
      UseOnePtAggregationAlgorithm = 0   [unused]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [unused]
      UseMaxLinkAggregationAlgorithm = 1   [unused]
      UseIsolatedNodeAggregationAlgorithm = 1   [unused]
      UseEmergencyAggregationAlgorithm = 1   [unused]
      aggregation: preserve Dirichlet points = 0   [unused]
      aggregation: enable phase 1 = 1   [unused]
      aggregation: enable phase 2a = 1   [unused]
      aggregation: enable phase 2b = 1   [unused]
      aggregation: enable phase 3 = 1   [unused]
      OnePt aggregate map name = 
      
      Build (MueLu::AmalgamationFactory)
      [empty list]
      
      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace
      
      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]
      
     [empty list]
     
    Damping factor = 1.33333
    
    Transpose P (MueLu::TransPFactory)
    [empty list]
    
   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]
   
   Keep AP Pattern = 0
   Keep RAP Pattern = 0
   implicit transpose = 0
   CheckMainDiagonal = 0
   RepairMainDiagonal = 0
   
  startLevel = 1
  minRowsPerProcessor = 2000
  nonzeroImbalance = 1.327
  remapPartitions = 1
  numRemapValues = 4
  alwaysKeepProc0 = 1
  repartition: print partition distribution = 0
  
 type = Interpolation
 implicit = 0
 implicit transpose = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 implicit = 0
 implicit transpose = 0
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 
 Setup Smoother (MueLu::AmesosSmoother{type = Superlu})
 presmoother -> 
  [empty list]
 
 
 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 4
 Operator complexity = 1.48
 Max Coarse Size     = 1000
 Implicit Transpose  = false
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  4
 A 1    3335  10015     3.00  1
 A 2    1111   3331     3.00  1
 A 3     371   1111     2.99  1
 
 Smoother (level 0) both : MueLu::IfpackSmoother{type = Chebyshev}
 
 Smoother (level 1) both : MueLu::IfpackSmoother{type = Chebyshev}
 
 Smoother (level 2) both : MueLu::IfpackSmoother{type = Chebyshev}
 
 Smoother (level 3) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 3) post : no smoother
 
