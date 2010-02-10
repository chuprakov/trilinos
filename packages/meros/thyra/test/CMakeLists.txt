INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_handles
  SOURCES test_handles.cpp
  ARGS --verbose
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_block_op
  SOURCES test_block_op.cpp
  ARGS --verbose
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_linear_combination
  SOURCES test_linear_combination
  ARGS --verbose
  NUM_MPI_PROCS 1
  )


ASSERT_DEFINED(CMAKE_HOST_SYSTEM_NAME)
IF (CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  SET(testEpetraOperatorWrapper_tol 1e-4)
ELSE()
  SET(testEpetraOperatorWrapper_tol 1e-3)
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  testEpetraOperatorWrapper
  SOURCES
    testEpetraOperatorWrapper.cpp
    ${COMMON_SOURCES}
  ARGS --tol=${testEpetraOperatorWrapper_tol}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
