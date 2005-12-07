#!/usr/bin/env sh
#
# This script performs the refactoring from the old MOOCHO to MOOCHO
# as a Trilinos package proper.  This script should be run from the base directory of
# any source tree that you want to change.
#
# To run this script, you must have the variable $MOOCHO_BASE_DIR set
# to the base MOOCHO directory so that
#
#    $MOOCHO_BASE_DIR/Moocho/src/MoochoUtilities/refactoring
#
# Trilinos/commonTools/refactoring must be added to your path.
#

$MOOCHO_BASE_DIR/moocho/src/MoochoUtilities/refactoring/from-old-to-trilinos.20051206.sh
