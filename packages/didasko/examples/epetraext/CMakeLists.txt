
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../epetraext/src/zoltan)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraEXT_ex1
  SOURCES ex1.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )
