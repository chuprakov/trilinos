
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  testInverseApply
  SOURCES testInverseApply.cpp
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
