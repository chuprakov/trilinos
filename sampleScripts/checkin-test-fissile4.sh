#!/bin/bash

# Used to test Trilinos on any of the ORNL fissle 4 machines
# (e.g. u233, u235, pu239, and pu241).

# NOTE: To use this, you must first prepend /opt/trilinos-toolset/bin
# to your path to find eg and cmake!

# NOTE: This script automatically picks up any CASL VRI related extra
# repos and adds them to --extra-repos.  If you want to override that,
# you can just pass in --extra-repos=??? to drop off extra repos or
# select the set that you want.

EXTRA_ARGS=$@

# The default location for this directory tree is:
#
#  Trilinos.base
#    Trilinos
#    BUILDS
#      CHECKIN
#
if [ "$TRILINOS_BASE_DIR" == "" ] ; then
  TRILINOS_BASE_DIR=../..
fi

TRILINOS_BASE_DIR_ABS=$(readlink -f $TRILINOS_BASE_DIR)

TRILINOS_TOOLSET_BASE=/opt/gcc-4.5.1/trilinos-toolset

EXTRA_REPOS_FULL_LIST="LIMEExt PSSDriversExt"

echo "
-DTrilinos_EXTRA_LINK_FLAGS:STRING='-Wl,-rpath,$TRILINOS_TOOLSET_BASE/lib64'
-DTPL_BLAS_LIBRARIES=/usr/lib64/libblas.so.3
-DTPL_LAPACK_LIBRARIES=/usr/lib64/liblapack.so.3
-DMPI_BASE_DIR:PATH=$TRILINOS_TOOLSET_BASE
" > MPI_DEBUG.config

echo "
-DTrilinos_EXTRA_LINK_FLAGS:STRING='-Wl,-rpath,$TRILINOS_TOOLSET_BASE/lib64'
-DTPL_BLAS_LIBRARIES=/usr/lib64/libblas.so.3
-DTPL_LAPACK_LIBRARIES=/usr/lib64/liblapack.so.3
-DCMAKE_CXX_COMPILER:PATH=$TRILINOS_TOOLSET_BASE/bin/g++
-DCMAKE_C_COMPILER:PATH=$TRILINOS_TOOLSET_BASE/bin/gcc
" > SERIAL_RELEASE.config

#
# Extra intel builds added with --extra-builds=INTEL_11064_SERIAL_DEBUG,...
#
# NOTE: You must do 'source /opt/casldev/env/casl_dev_env.sh' before
# using the intel builds.

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${TRILINOS_BASE_DIR_ABS}/Trilinos/cmake/ctest/drivers/pu241/intel-11.064-options.cmake
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON
-DTrilinos_ENABLE_TESTS:BOOL=ON
-DDART_TESTING_TIMEOUT:STRING=180.0
" > INTEL_11064_SERIAL_DEBUG.config


#
# Load up the list of extra repos based on what is present:
#

EXTRA_REPOS=
for extra_repo in $EXTRA_REPOS_FULL_LIST; do
  #echo $extra_repo
  EXTRA_REPO_PATH=$TRILINOS_BASE_DIR/Trilinos/$extra_repo
  #echo $EXTRA_REPO_PATH
  if [ -d $EXTRA_REPO_PATH ]; then
    EXTRA_REPOS=$EXTRA_REPOS$extra_repo,
  fi
done
#echo "EXTRA_REPOS=$EXTRA_REPOS"

#
# Invocation
#

$TRILINOS_BASE_DIR/Trilinos/checkin-test.py \
--extra-repos=$EXTRA_REPOS \
-j16 \
--ctest-timeout=180 \
$EXTRA_ARGS  

# NOTE: By default we use 16 processes which is 1/2 of the 32
# processes on this machine.  This way two people can build and test
# Trilinos without taxing the machine too much.
