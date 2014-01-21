repartition: enable = 1
repartition: remap parts = 0
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
repartition: partitioner = zoltan   [default]

Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 relaxation: type = Symmetric Gauss-Seidel
 relaxation: sweeps = 1
 relaxation: damping factor = 1
 relaxation: zero starting solution = 1   [default]
 relaxation: backward mode = 0   [default]
 relaxation: use l1 = 0   [default]
 relaxation: l1 eta = 1.5   [default]
 relaxation: min diagonal value = 0   [default]
 relaxation: fix tiny diagonal entries = 0   [default]
 relaxation: check diagonal entries = 0   [default]
 relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
       lightweight wrap = 1
       aggregation threshold = 0
       Dirichlet detection threshold = 0
       algorithm = original
       
      Ordering = 0   [default]
      MaxNeighAlreadySelected = 0   [default]
      MinNodesPerAggregate = 2   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
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
      Fine level nullspace = Nullspace
      
      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      
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
   
  startLevel = 2
  minRowsPerProcessor = 800
  nonzeroImbalance = 1.2
  remapPartitions = 0
  numRemapValues = 4   [unused]
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
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 relaxation: type = Symmetric Gauss-Seidel
 relaxation: sweeps = 1
 relaxation: damping factor = 1
 relaxation: zero starting solution = 1   [default]
 relaxation: backward mode = 0   [default]
 relaxation: use l1 = 0   [default]
 relaxation: l1 eta = 1.5   [default]
 relaxation: min diagonal value = 0   [default]
 relaxation: fix tiny diagonal entries = 0   [default]
 relaxation: check diagonal entries = 0   [default]
 relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
       lightweight wrap = 1
       aggregation threshold = 0
       Dirichlet detection threshold = 0
       algorithm = original
       
      Ordering = 0   [default]
      MaxNeighAlreadySelected = 0   [default]
      MinNodesPerAggregate = 2   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
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
      Fine level nullspace = Nullspace
      
      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      
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
   
  startLevel = 2
  minRowsPerProcessor = 800
  nonzeroImbalance = 1.2
  remapPartitions = 0
  numRemapValues = 4   [unused]
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
 
 Setup Smoother (MueLu::Amesos2Smoother{type = Superlu})
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
 
 Smoother (level 0) both : "Ifpack2::Relaxation": { MatrixType: "Tpetra::CrsMatrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > >", Status: initialized, computed, "relaxation: type": Symmetric Gauss-Seidel, "relaxation: sweeps": 1, "relaxation: damping factor": 1, "Global number of rows": 9999, "Global number of columns": 9999 }
 
 Smoother (level 1) both : "Ifpack2::Relaxation": { MatrixType: "Tpetra::CrsMatrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > >", Status: initialized, computed, "relaxation: type": Symmetric Gauss-Seidel, "relaxation: sweeps": 1, "relaxation: damping factor": 1, "Global number of rows": 3333, "Global number of columns": 3333 }
 
 Smoother (level 2) pre  : SuperLU solver interface
 Smoother (level 2) post : no smoother
 
