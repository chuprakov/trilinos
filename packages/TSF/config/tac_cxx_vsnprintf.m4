dnl @synopsis AC_CXX_VSNPRINTF
dnl
dnl If the compiler recognizes vsnprintf as a function for IO,
dnl define HAVE_VSNPRINTF.  If this test fails, use KL's workaround
dnl for vsprintf.
dnl Note that we try to compile two versions of this routine, one using cstdio and
dnl another using stdio.h.  This approach is used to eliminate the need to test which
dnl of the two header files is present.  If one or both is usable the test will return true.
dnl

AC_DEFUN([AC_CXX_VSNPRINTF],
[AC_CACHE_CHECK([[whether the compiler recognizes vsnprintf as supported IO function]],
ac_cv_cxx_vsnprintf,
[ AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <cstdio>
#include <cstdarg>
#include <string>

void myTest(const char* f, ...)
{
  va_list args;  
	va_start(args, f); 

  char str[100];
  vsnprintf(str, 99, f, args);
}
],
[
  myTest("my test %d", 1);
],
 ac_cv_cxx_vsnprintf1=yes, ac_cv_cxx_vsnprintf1=no)

AC_TRY_COMPILE([
#include <stdio.h>
#include <stdarg.h>
#include <string>


void myTest(const char* f, ...)
{
  va_list args;  
	va_start(args, f); 

  char str[100];
  vsnprintf(str, 99, f, args);
}
],
[
  myTest("my test %d", 1);
],
 ac_cv_cxx_vsnprintf2=yes, ac_cv_cxx_vsnprintf2=no)

if (test "$ac_cv_cxx_vsnprintf1" = yes || test "$ac_cv_cxx_vsnprintf2" = yes); then
 ac_cv_cxx_vsnprintf=yes
else
 ac_cv_cxx_vsnprintf=no
fi
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_vsnprintf" = yes; then
  AC_DEFINE(HAVE_VSNPRINTF,,[define if vsnprintf is supported])
fi
])
