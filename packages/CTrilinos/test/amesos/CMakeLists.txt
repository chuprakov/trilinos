INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)


INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR}/.. )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CAmesos_UnitTests
  SOURCES
    CAmesos_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  CAmesos_BaseSolver_UnitTests
#  SOURCES
#    CAmesos_BaseSolver_UnitTests.cpp
#    ../CTrilinos_UnitTestHelpers.cpp
#    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#  )

