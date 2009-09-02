
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Phdmesh)

IF(NOT HAVE_SYS_TIME_H)
  SET(PHDMESH_NO_SYS_TIME TRUE)
ELSE()
  SET(PHDMESH_NO_SYS_TIME FALSE)
ENDIF()

#
# B) Set up package-specific options
#

IF(TPL_ENABLE_ExodusII)
  SET(${PACKAGE_NAME}_ENABLE_EXODUSII_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_EXODUSII_DEFAULT OFF)
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EXODUSII
  PHDMESH_HAS_SNL_EXODUSII
  "Enable ExodusII support for phdmesh mesh IO."
  ${Phdmesh_ENABLE_EXODUSII_DEFAULT} )

#convert the mpi option to the one that is used in the
#config header file
IF(TPL_ENABLE_MPI)
  SET(PHDMESH_HAS_MPI TRUE)
ENDIF()


#Checking for the pthread library. Not sure why there is
#a positive and negative define though.
IF(${PACKAGE_NAME}_ENABLE_Pthread)
  MESSAGE("--  Building with thread support")
  SET(HAVE_PTHREAD TRUE)
ELSE()
  MESSAGE("--  Building without thread support")
  SET(NO_PTHREADS TRUE)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(tests)


#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
