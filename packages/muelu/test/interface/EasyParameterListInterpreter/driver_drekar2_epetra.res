verbosity = test
number of equations = 1
coarse: max size = 1000
sa: use filtered matrix = 1
filtered matrix: use lumping = 1
smoother: type = CHEBYSHEV
aggregation: drop scheme = laplacian
aggregation: drop tol = 0.02
repartition: enable = 1
repartition: min rows per proc = 2000
repartition: max imbalance = 1.327
repartition: start level = 1
repartition: remap parts = 1
repartition: keep proc 0 = 1
repartition: partitioner = zoltan2
max levels = 10   [default]
debug: graph level = -1   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
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
 chebyshev: max eigenvalue = -1   [default]
 chebyshev: min diagonal value = 0   [default]
 chebyshev: operator inv diagonal = 0   [default]
 chebyshev: use block mode = 0   [default]
 chebyshev: solve normal equations = 0   [default]

Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      Dirichlet detection threshold = 0
      aggregation threshold = 0.02
      algorithm = laplacian

     lumping = 1

     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      Ordering = 0   [default]
      MaxNeighAlreadySelected = 0   [default]
      MinNodesPerAggregate = 2   [default]
      UseOnePtAggregationAlgorithm = 1   [default]
      UseSmallAggregatesAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]
      SmallAgg aggregate map name =    [default]

      Build (MueLu::AmalgamationFactory)
      [empty list]

      Nullspace factory (MueLu::NullspaceFactory)
      [empty list]

      Build (MueLu::CoarseMapFactory)
      [empty list]

     [empty list]

    Damping factor = 1.33333

    Transpose P (MueLu::TransPFactory)
    [empty list]

   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]

   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]

  startLevel = 1
  minRowsPerProcessor = 2000
  nonzeroImbalance = 1.327
  remapPartitions = 1
  alwaysKeepProc0 = 1

 type = Interpolation
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 chebyshev: max eigenvalue = -1   [default]
 chebyshev: min diagonal value = 0   [default]
 chebyshev: operator inv diagonal = 0   [default]
 chebyshev: use block mode = 0   [default]
 chebyshev: solve normal equations = 0   [default]

Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      Dirichlet detection threshold = 0
      aggregation threshold = 0.02
      algorithm = laplacian

     lumping = 1

     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      Ordering = 0   [default]
      MaxNeighAlreadySelected = 0   [default]
      MinNodesPerAggregate = 2   [default]
      UseOnePtAggregationAlgorithm = 1   [default]
      UseSmallAggregatesAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]
      SmallAgg aggregate map name =    [default]

      Build (MueLu::AmalgamationFactory)
      [empty list]

      Nullspace factory (MueLu::NullspaceFactory)
      [empty list]

      Build (MueLu::CoarseMapFactory)
      [empty list]

     [empty list]

    Damping factor = 1.33333

    Transpose P (MueLu::TransPFactory)
    [empty list]

   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]

   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]

  startLevel = 1
  minRowsPerProcessor = 2000
  nonzeroImbalance = 1.327
  remapPartitions = 1
  alwaysKeepProc0 = 1

 type = Interpolation
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 chebyshev: max eigenvalue = -1   [default]
 chebyshev: min diagonal value = 0   [default]
 chebyshev: operator inv diagonal = 0   [default]
 chebyshev: use block mode = 0   [default]
 chebyshev: solve normal equations = 0   [default]

Level 3
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      Dirichlet detection threshold = 0
      aggregation threshold = 0.02
      algorithm = laplacian

     lumping = 1

     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      Ordering = 0   [default]
      MaxNeighAlreadySelected = 0   [default]
      MinNodesPerAggregate = 2   [default]
      UseOnePtAggregationAlgorithm = 1   [default]
      UseSmallAggregatesAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]
      SmallAgg aggregate map name =    [default]

      Build (MueLu::AmalgamationFactory)
      [empty list]

      Nullspace factory (MueLu::NullspaceFactory)
      [empty list]

      Build (MueLu::CoarseMapFactory)
      [empty list]

     [empty list]

    Damping factor = 1.33333

    Transpose P (MueLu::TransPFactory)
    [empty list]

   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]

   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]

  startLevel = 1
  minRowsPerProcessor = 2000
  nonzeroImbalance = 1.327
  remapPartitions = 1
  alwaysKeepProc0 = 1

 type = Interpolation
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
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
 A 0    9999  29995     3.00  1
 A 1    3333   9997     3.00  1
 A 2    1111   3331     3.00  1
 A 3     371   1111     2.99  1

 Smoother (level 0) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 1) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 2) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 3) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 3) post : no smoother

