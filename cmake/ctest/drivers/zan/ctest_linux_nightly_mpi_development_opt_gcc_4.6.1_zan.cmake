
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.zan.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV2)
SET(COMPILER_VERSION "GCC-4.6.1")
SET(ENV{LD_LIBRARY_PATH} "/home/jmwille/install/mpc-0.9/lib:/home/jmwille/install/mpfr-2.4.2/lib:/home/jmwille/install/gmp-4.3.2/lib:/home/jmwille/install/gcc-4.6.1/lib64:$ENV{LD_LIBRARY_PATH}")
SET(CTEST_TEST_TYPE EXPERIMENTAL)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
#Stokhos is explicitly disabled below to prevent the package from being
#implicitly enabled.  Sundance depends on Stokhos.
#SET(EXTRA_EXCLUDE_PACKAGES Phalanx Stokhos Sundance)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON"
  "-DMPI_BASE_DIR:PATH=/home/jmwille/install/gcc-4.6.1/openmpi-1.4.3"
  "-DBoost_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.4.4/boost_1_46_1"
  "-DNetcdf_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.4.4/netcdf-4.1.3/include"
  "-DNetcdf_LIBRARY_DIRS=/home/trilinos/tpl/gcc4.4.4/netcdf-4.1.3/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
