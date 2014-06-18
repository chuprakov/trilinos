max levels = 1
coarse: type = RELAXATION
verbosity = test
coarse: max size = 2000   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
aggregation: visualize = 0   [default]
repartition: enable = 0   [default]

Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 relaxation: type = Jacobi   [default]
 relaxation: sweeps = 1   [default]
 relaxation: damping factor = 1   [default]
 relaxation: zero starting solution = 1   [default]
 relaxation: backward mode = 0   [default]
 relaxation: use l1 = 0   [default]
 relaxation: l1 eta = 1.5   [default]
 relaxation: min diagonal value = 0   [default]
 relaxation: fix tiny diagonal entries = 0   [default]
 relaxation: check diagonal entries = 0   [default]
 relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
 
 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 1
 Operator complexity = 1.00
 Max Coarse Size     = 2000
 Implicit Transpose  = false
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1
 
 Smoother (level 0) pre  : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Jacobi, sweeps: 1, damping factor: 1, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
 Smoother (level 0) post : no smoother
 
