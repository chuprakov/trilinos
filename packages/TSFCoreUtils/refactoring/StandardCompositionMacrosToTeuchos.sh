#!/bin/sh
# Warning, to use this script you must first define TRILINOS_HOME and add the
# base path for the scrpts 'token-replace-list-r' and 'string-replace-list-r' to your path.
# This script can be run over and over again on the same files without harm.
string-replace-list-r $TRILINOS_HOME/packages/TSFCoreUtils/refactoring/StandardCompositionMacrosToTeuchos-string-replace.list

