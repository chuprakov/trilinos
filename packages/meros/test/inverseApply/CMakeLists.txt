
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  estInverseApply
  SOURCES testInverseApply.cpp
  ARGS 
  COMM serial mpi
  )
