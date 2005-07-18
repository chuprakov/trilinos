#!/bin/sh
#
# This script can be used to automatically upgrade code for the refactoring of code
# from TSFCore to Thyra.  This is for packages that are not in the TSFCore namespace.
# To use this script you must set the environment
# variable TRILINOS_HOME to the path of your Trilinos base directory and you
# must add $TRILINO_HOME/commonTools/refactoring to your path.

# This script only changes the names of classes and header files and does not
# make all of the changes that are needed.  A few changes that are not automated are:
#
# 1) Change A.range()->createMember() to createMember(A.range()) and do the same for
#    all calls of createMember() and createMembers().  This was needed to insure
#    that vectors and multi-vectors live past the vector space that creates them.
#
# 2) Change A.apply(NOTRANS,...) to Thyra::apply(A,NOTRANS,...) and A.opSupported(NOTRANS) to
#    Thyra::opSupported(A,NOTRANS).  This was needed to preserve use by code that does not
#    use more than one scalar type but the new Thyra::LinearOpBase interface is not
#    templated on RangeScalar and DomainScalar types.
#
# 3) Most existing TSFCore::LinearOp subclasses should inherit from Thyra::SingleScalarLinearOpBase
#    if they implemented multi-vector version of apply(...) or from Thyra::SingleRhsLinearOpBase
#    if they implemented single-vector version of apply(...).
#
# 4) Code that directly used TSFCore::EpetraVector or TSFCore::EpetraVectorSpace now needs to
#    do something else or use the new wrapper functions in Thyra_EpetraThyraWrappers.hpp.
#
# More changes may be needed other than these.
#


# Run this script from the base directory of any code that you would like to
# upgrade.  This scrpit can be run mutiple times with no side effects.

if [ -z "$TRILINOS_HOME" ]; then 
  echo "$0 Usage:"
  echo "Please set TRILINOS_HOME env variable"
  echo "Please add $TRILINOS_HOME/commonTools/refactoring to your path"
  exit 1
fi

token-replace-list-r $TRILINOS_HOME/packages/TSFCore/refactoring/TSFCore-to-Thyra.out-TSFCore.token-replace.list
string-replace-list-r $TRILINOS_HOME/packages/TSFCore/refactoring/TSFCore-to-Thyra.string-replace.list
