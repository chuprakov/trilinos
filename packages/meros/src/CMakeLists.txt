
INCLUDE(PackageLibraryMacros)
INCLUDE(CheckIncludeFileCXX)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

#
# Core files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Meros_AztecSolveStrategy.hpp
  Meros_ConfigDefs.h
  Meros_IdentityOperator.hpp
  Meros_InverseOperator.hpp
  Meros_LinearSolver.hpp
  Meros_LinearSolveStrategy.hpp
  Meros_LSCOperatorSource.h
  Meros_LSCPreconditionerFactory.h
  Meros_PCDOperatorSource.h
  Meros_PCDPreconditionerFactory.h
  Meros_Preconditioner.hpp
  Meros_PreconditioningStrategy.hpp
  Meros_SIMPLEOperatorSource.h
  Meros_SIMPLEPreconditionerFactory.h
  Meros_ZeroOperator.hpp
  )

APPEND_SET(SOURCES
  Meros_LSCOperatorSource.cpp
  Meros_LSCPreconditionerFactory.cpp
  Meros_PCDOperatorSource.cpp
  Meros_PCDPreconditionerFactory.cpp
  Meros_SIMPLEOperatorSource.cpp
  Meros_SIMPLEPreconditionerFactory.cpp
  )

#
# util headers
#

APPEND_SET(HEADERS
  Meros_Version.h
  )


######################################

APPEND_SET(HEADERS
  )

APPEND_SET(SOURCES
  )

######################################
IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  meros
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
