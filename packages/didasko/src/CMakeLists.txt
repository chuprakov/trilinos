
#This file is only here to set up the config file in the right spot.
#potentially we could put the latex building logic in here as well at
#some point.

INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PRINT_VAR(HAVE_CASSERT)
PRINT_VAR(CMATH)
PRINT_VAR(CSTDIO)
PRINT_VAR(HAVE_CSTDLIB)
PRINT_VAR(HAVE_IOSTREAM)
PRINT_VAR(STRING)
PRINT_VAR(UNISTD_H)
PRINT_VAR(HAVE_ZOLTAN)

IF(HAVE_DIDASKO_EPETRAEXT AND Trilinos_ENABLE_Zoltan)
  SET(HAVE_ZOLTAN ON)
ENDIF()

PRINT_VAR(HAVE_ZOLTAN)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

