dnl @synopsis TAC_ARG_CONFIG_BOUNDSCHECK
dnl
dnl Support for boundschecking in TSFArray
dnl
dnl --enable-tsf-boundscheck       - Turns on boundschecking
dnl
dnl @author Kevin Long (krlong@sandia.gov)
dnl
AC_DEFUN([TAC_ARG_CONFIG_BOUNDSCHECK],
[
AC_ARG_ENABLE(boundschecking,
AC_HELP_STRING([--enable-boundschecking],[enable TSFArray boundschecking]),
[
BNDSCHK=yes
],
[BNDSCHK=no]
)

AC_MSG_CHECKING(whether we are using TSF array boundschecking)
AC_MSG_RESULT([${BNDSCHK}])

if test "X$BNDSCHK" = "Xyes"; then
  AC_DEFINE(HAVE_ARRAY_BOUNDSCHECK,,[Define if we are using array boundschecking])
fi

])

