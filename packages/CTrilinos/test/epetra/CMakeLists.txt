INCLUDE(PackageAddExecutableAndTest)


#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CEpetra_UnitTestsHelpers_UnitTests
  SOURCES
    CEpetra_UnitTestHelpers_UnitTests.cpp
    CEpetra_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Comm_UnitTests
  SOURCES
    CEpetra_Comm_UnitTests.cpp
    CEpetra_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CEpetra_SerialComm_UnitTests
  SOURCES
    CEpetra_SerialComm_UnitTests.cpp
    CEpetra_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  CEpetra_MpiComm_UnitTests
#  SOURCES
#    CEpetra_MpiComm_UnitTests.cpp
#    CEpetra_UnitTestHelpers.cpp
#    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
#  COMM mpi
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Map_UnitTests
  SOURCES
    CEpetra_Map_UnitTests.cpp
    CEpetra_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Vector_UnitTests
  SOURCES
    CEpetra_Vector_UnitTests.cpp
    CEpetra_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

