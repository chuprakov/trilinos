dnl @synopsis TAC_ARG_CHECK_FEATURE_SUB(FEATURE_NAME, SUB_FEATURE_NAME, HAVE_NAME)
dnl
dnl Test for --enable-${FEATURE_NAME} and set to yes
dnl (unless --disable-default-packages is specified) if not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help define whether code dependent upon packages should be built by 
dnl default (unless --disable-default-packages is specified) should be built.
dnl  For example:
dnl
dnl TAC_ARG_CHECK_FEATURE_SUB(tsfcore, epetra, TSFCORE_EPETRA)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined and --disable-default-packages is not 
dnl specified then HAVE_EPETRA will be defined, if --enable-epetra is not 
dnl defined to be "yes" and --disable-default-packages is specified, 
dnl HAVE_EPETRA will not be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on Trilinos/config/tac_arg_enable_default_feature.m4
dnl @author Heidi Thornquist <hkthorn@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CHECK_FEATURE_SUB],
[
AC_ARG_ENABLE([$1-$2],, ac_cv_use_$1_$2=$enableval, ac_cv_use_$1_$2=no)

if test "X$ac_cv_use_$1_$2" != "Xno"; then
  AC_DEFINE([HAVE_$3],,[Define if want to build $1-$2])
fi
])

