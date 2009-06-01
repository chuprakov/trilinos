
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simpleEpetraExternCInterface_test  
  SOURCES simpleEpetraExternCInterface.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

SET_SOURCE_FILES_PROPERTIES( simpleEpetraExternCInterface.c PROPERTIES LANGUAGE CXX )
