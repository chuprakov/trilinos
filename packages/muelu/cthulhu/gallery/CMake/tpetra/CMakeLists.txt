INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(MatrixDriver_SOURCES
  MueLu_MatrixDriver.cpp
  MueLu_MatrixFactory.hpp
  MueLu_MatrixTypes.hpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixDriver_Tpetra
  DIRECTORY ../../
  SOURCES ${MatrixDriver_SOURCES}
  DEFINES -DCTHULHU_USE_TPETRA
  )
