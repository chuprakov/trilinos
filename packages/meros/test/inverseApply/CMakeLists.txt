
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  testInverseApply
  SOURCES testInverseApply.cpp
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
