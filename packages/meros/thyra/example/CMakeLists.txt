INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/../thyra/example/operator_vector)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_mpi
  SOURCES sillyCgSolve_mpi.cpp
  )
