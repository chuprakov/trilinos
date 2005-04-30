#!/bin/sh
#

# This script can be used to automatically upgrade code for the move of code
# in TSFCoreUtils into Thyra.  To use this script you must set the environment
# variable TRILINOS_HOME to the path of your Trilinos base directory and you
# must add $TRILINO_HOME/commonTools/refactoring to your path.

# Run this script from the base directory of any code that you would like to
# upgrade.  This scrpit can be run mutiple times with no side effects.

token-replace-list-r $TRILINOS_HOME/packages/TSFCoreUtils/refactoring/TSFCoreUtils-to-Thyra.token-replace.list
string-replace-list-r $TRILINOS_HOME/packages/TSFCoreUtils/refactoring/TSFCoreUtils-to-Thyra.string-replace.list
