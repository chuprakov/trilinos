
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  saddle_lsc
  SOURCES saddle_lsc.cpp
  ARGS 
  COMM serial mpi
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  saddle_lsc_epetra
  SOURCES saddle_lsc_epetra.cpp
  ARGS 
  COMM serial mpi
  )
