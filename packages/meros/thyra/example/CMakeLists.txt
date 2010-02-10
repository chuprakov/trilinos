INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/../thyra/example/operator_vector)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_mpi
  SOURCES sillyCgSolve_mpi.cpp
  )
