INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST( ExampleNLPBanded
  SOURCES
    ExampleNLPBandedMain.cpp
    NLPInterfacePack_ExampleNLPBanded.cpp
  ARGS
    "--nD=2 --nI=1 --bw=1 --mI=0 --xo=0.1 --diag-scal=1e+1 --nlp-selects-basis"
    "--nD=2 --nI=2 --bw=2 --mI=0 --xo=0.1 --diag-scal=1e+1"
    "--nD=1000 --nI=8 --bw=5 --mI=0 --nlp-selects-basis --xo=0.1 --diag-scal=1e+3"
    "--nD=1000 --nI=300 --bw=20 --mI=0 --nlp-selects-basis --xo=0.1 --diag-scal=1e+2"
    "--nD=8 --nI=2 --bw=2 --mI=4 --xo=0.1 --diag-scal=1e+1"
    "--nD=100 --nI=50 --bw=100 --mI=100 --nlp-selects-basis --xo=0.1 --diag-scal=1e+1"
    "--nD=100 --nI=4 --bw=20 --mI=0 --xo=0.1 --diag-scal=1e+3"
    "--nD=4 --nI=4 --bw=2 --mI=0 --xo=0.1 --diag-scal=1e+3 --xIl=0.01"
    "--nD=2 --nI=2 --ignore-constraints --xo=0.1 --xIl=0.01"
    "--nD=100 --nI=100 --ignore-constraints --xo=0.1 --xIl=0.01"
  NUM_MPI_PROCS 1
  )


