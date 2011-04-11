#!/bin/bash

echo
echo "Starting nightly Trilinos testing on pu241: `date`"
echo

setenv PATH "${PATH}:/opt/casldev/env"
eval `vera_dev_env.py load casl`
eval `vera_dev_env.py load intel/12`

BASEDIR=/home/casl-vri-admin/Dashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/pu241
TRILINOS_REPOSITORY_LOCATION="cgbaker@software.sandia.gov:/space/git/Trilinos"

export TDD_PARALLEL_LEVEL=2
export TDD_CTEST_TEST_TYPE=Experimental

export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/CDash/submit.php?project=TrilinosDriver"

time env python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on pu241: `date`"
echo

echo "Finished nightly Trilinos CMake tests pu241: http://casl-dev.ornl.gov/CDash/index.php?project=Trilinos" | mailx -s "Nightly CTest: pu241" casl-vri-admin@localhost
