  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for trilinos-test using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET_DEFAULT( CTEST_BUILD_FLAGS "-j7 -i" )

  SET_DEFAULT( CTEST_PARALLEL_LEVEL "8" )

#  SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
  SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

  # Only turn on PyTrilinos for shared libraries
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} PyTrilinos TriKota Optika)
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
    "-DMEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
    "-DMesquite_ENABLE_TESTS:BOOL=ON"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DNetcdf_LIBRARY_DIRS=/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/lib"
    "-DNetcdf_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/include"
    "-DTPL_ENABLE_MATLAB=OFF"
    )

  IF (BUILD_TYPE STREQUAL "DEBUG")
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
      )
  ENDIF()

  SET_DEFAULT(COMPILER_VERSION "PGI-11.1")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      )
#      "-DMPI_BASE_DIR:PATH=/home/trilinos/openmpi-1.4"
#      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_trilinos-test_openmpi_1.2.7.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DCMAKE_CXX_COMPILER:FILEPATH=/home/trilinos/pgi11.1/linux86-64/11.1/bin/pgCC"
      "-DCMAKE_C_COMPILER:FILEPATH=/home/trilinos/pgi11.1/linux86-64/11.1/bin/pgcc"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=/home/trilinos/pgi11.1/linux86-64/11.1/bin/pgf90"
      "-DTrilinos_EXTRA_LINK_FLAGS:STRING='-L/home/trilinos/pgi11.1/linux86-64/11.1/lib -L/usr/lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl  -lnspgc -lpgc  -lrt -lpthread  -lm -lgcc -lc -lgcc'"
      "-DTrilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON"
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_trilinos-test_gcc-4.1.2.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
