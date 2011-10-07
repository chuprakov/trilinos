#
# This file is meant to be included in the configuration of Trilinos to
# disable packages and other critical varaibles internally.  These excludes
# can be overwritten on the command line in the CMake cache.
#
# These options are set for *all* compilers for all builds so they have to be
# 100% general!
#

# Must include first so that it defines Trilinos_EXCLUDE_PACKAGES
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-exclude-trilinos-packages.cmake)
#PRINT_VAR(Trilinos_EXCLUDE_PACKAGES)

# Put in hard disables 
FOREACH(EXCLUDED_PACKAGE ${Trilinos_EXCLUDE_PACKAGES})
  SET(${PROJECT_NAME}_ENABLE_${EXCLUDED_PACKAGE} OFF CACHE BOOL
    "Disabled in casl-core-enables-disables.cmake")
  #SET(${PROJECT_NAME}_ENABLE_${EXCLUDED_PACKAGE} OFF)
  #PRINT_VAR(${PROJECT_NAME}_ENABLE_${EXCLUDED_PACKAGE})
ENDFOREACH()

# Turn off float and complex testing because CASL does not need them
SET(Teuchos_ENABLE_FLOAT OFF CACHE BOOL "")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL "")

# Turn off STK tests since they are constantly failing.  NOTE: Since CASL is
# not developing on STK, only using it, this should not represent a big risk
# for STK or CASL.
SET(STK_ENABLE_TESTS OFF CACHE BOOL "")
SET(STK_ENABLE_EXAMPLES OFF CACHE BOOL "")

# Turn on configure timing
SET(Trilinos_ENABLE_CONFIGURE_TIMING ON CACHE BOOL "")

# Don't create *Config.cmake files since they are massively expensive to create
SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL "")
