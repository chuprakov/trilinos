
#This file is only here to set up the config file in the right spot.
#potentially we could put the latex building logic in here as well at
#some point.

INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

IF(HAVE_DIDASKO_EPETRAEXT AND Trilinos_ENABLE_Zoltan)
  SET(HAVE_ZOLTAN ON)
ENDIF()

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

