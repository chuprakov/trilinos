verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
smoother: pre or post = both   [default]
smoother: type = RELAXATION   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
repartition: enable = 0   [default]
smoother: pre params -> 
 relaxation: sweeps = 3   [unused]
smoother: post params -> 
 relaxation: sweeps = 0   [unused]

Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 presmoother -> 
  relaxation: sweeps = 3
  relaxation: type = Jacobi   [default]
  relaxation: damping factor = 1   [default]
  relaxation: zero starting solution = 1   [default]
  relaxation: backward mode = 0   [default]
  relaxation: use l1 = 0   [default]
  relaxation: l1 eta = 1.5   [default]
  relaxation: min diagonal value = 0   [default]
  relaxation: fix tiny diagonal entries = 0   [default]
  relaxation: check diagonal entries = 0   [default]
  relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 postsmoother -> 
  relaxation: sweeps = 0
  relaxation: type = Jacobi   [default]
  relaxation: damping factor = 1   [default]
  relaxation: zero starting solution = 1   [default]
  relaxation: backward mode = 0   [default]
  relaxation: use l1 = 0   [default]
  relaxation: l1 eta = 1.5   [default]
  relaxation: min diagonal value = 0   [default]
  relaxation: fix tiny diagonal entries = 0   [default]
  relaxation: check diagonal entries = 0   [default]
  relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
Level 1
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    Dirichlet detection threshold = 0
    aggregation threshold = 0
    algorithm = original
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
   UseOnePtAggregationAlgorithm = 1   [default]
   UseSmallAggregatesAggregationAlgorithm = 0   [default]
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
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 presmoother -> 
  relaxation: sweeps = 3
  relaxation: type = Jacobi   [default]
  relaxation: damping factor = 1   [default]
  relaxation: zero starting solution = 1   [default]
  relaxation: backward mode = 0   [default]
  relaxation: use l1 = 0   [default]
  relaxation: l1 eta = 1.5   [default]
  relaxation: min diagonal value = 0   [default]
  relaxation: fix tiny diagonal entries = 0   [default]
  relaxation: check diagonal entries = 0   [default]
  relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 postsmoother -> 
  relaxation: sweeps = 0
  relaxation: type = Jacobi   [default]
  relaxation: damping factor = 1   [default]
  relaxation: zero starting solution = 1   [default]
  relaxation: backward mode = 0   [default]
  relaxation: use l1 = 0   [default]
  relaxation: l1 eta = 1.5   [default]
  relaxation: min diagonal value = 0   [default]
  relaxation: fix tiny diagonal entries = 0   [default]
  relaxation: check diagonal entries = 0   [default]
  relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    Dirichlet detection threshold = 0
    aggregation threshold = 0
    algorithm = original
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
   UseOnePtAggregationAlgorithm = 1   [default]
   UseSmallAggregatesAggregationAlgorithm = 0   [default]
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
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 
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
 
 Smoother (level 0) pre  : "Ifpack2::Relaxation": { MatrixType: "Tpetra::CrsMatrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > >", Status: initialized, computed, "relaxation: type": Jacobi, "relaxation: sweeps": 3, "relaxation: damping factor": 1, "Global number of rows": 9999, "Global number of columns": 9999 }
 Smoother (level 0) post : "Ifpack2::Relaxation": { MatrixType: "Tpetra::CrsMatrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > >", Status: initialized, computed, "relaxation: type": Jacobi, "relaxation: sweeps": 0, "relaxation: damping factor": 1, "Global number of rows": 9999, "Global number of columns": 9999 }
 
 Smoother (level 1) pre  : "Ifpack2::Relaxation": { MatrixType: "Tpetra::CrsMatrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > >", Status: initialized, computed, "relaxation: type": Jacobi, "relaxation: sweeps": 3, "relaxation: damping factor": 1, "Global number of rows": 3333, "Global number of columns": 3333 }
 Smoother (level 1) post : "Ifpack2::Relaxation": { MatrixType: "Tpetra::CrsMatrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > >", Status: initialized, computed, "relaxation: type": Jacobi, "relaxation: sweeps": 0, "relaxation: damping factor": 1, "Global number of rows": 3333, "Global number of columns": 3333 }
 
 Smoother (level 2) pre  : SuperLU solver interface
 Smoother (level 2) post : no smoother
 
