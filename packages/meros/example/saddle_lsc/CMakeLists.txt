
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE(
  saddle_lsc
  SOURCES saddle_lsc.cpp
  COMM serial mpi
  )

TRILINOS_ADD_EXECUTABLE(
  saddle_lsc_epetra
  SOURCES saddle_lsc_epetra.cpp
  COMM serial mpi
  )

# 2008/12/20: rabartl: ToDo: Run these examples in order to test them.
# You can't just run them because these segfault when run without
# arguments!  Also, these are not being run from meros/test/definition
# so we should not be running them either.
