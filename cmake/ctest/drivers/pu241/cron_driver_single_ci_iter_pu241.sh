#!/bin/bash

echo
echo "Starting continuous integration Trilinos testing iteration on pu241: `date`"
echo

export PATH="${PATH}:/opt/casldev/env"
eval `vera_dev_env.py load casl`
eval `vera_dev_env.py load intel/12`
umask u=rwx,g=rwx,o=

BASEDIR=/home/casl-vri-admin/Dashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/pu241
TRILINOS_REPOSITORY_LOCATION="cgbaker@software.sandia.gov:/space/git/Trilinos"

export TDD_PARALLEL_LEVEL=4
export TDD_CTEST_TEST_TYPE=Continuous

# Submit the outer TDD tests to casl-dev always since these are CASL machines
export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/CDash/submit.php?project=TrilinosDriver"

time env python ../cron_driver.py

echo
echo "Ending continuous integration Trilinos testing iteration on pu241: `date`"
echo
