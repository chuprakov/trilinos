# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(TSFCoreUtils, 0.0, rabartl@sandia.gov)

# Hello World!
#  #np# - package name below
echo "----------------------------------------"
echo "Running TSFCoreUtils Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)
#  #np# - Change file names in next line 
# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

AM_CONFIG_HEADER(src/TSFCoreUtilsConfig.h:src/TSFCoreUtilsConfig.h.in)

# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
                                                                                
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)
                                                                                
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

#  #np# - Package name and version below
AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

#  #np# - Customize the checks for headers needed in package being added.
#NOTE: This section still needs to be updated for new_package, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([cmath math.h iomanip iomanip.h ostream iostream iostream.h],,,)
AC_CHECK_HEADERS([string string.h malloc.h stdexcept map map.h vector vector.h],,,)
AC_CHECK_HEADERS([memory memory.h typeinfo limits.h new iosfwd],,,)
AC_CHECK_HEADERS([sstream sstream.h strstream strstream.h],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  TSFCoreUtils requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  TSFCoreUtils requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  TSFCoreUtils requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  TSFCoreUtils requires at least one or the other.]])
fi

if test "X$ac_cv_header_ostream" != "Xyes" && test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  TSFCoreUtils requires at least one or the other.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  TSFCoreUtils requires at least one or the other.]])
fi
  
if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  TSFCoreUtils requires at least one or the other.]])
fi  

if test "X$ac_cv_header_malloc_h" != "Xyes"; then
 AC_FUNC_MALLOC
 #if test "X$ac_cv_func_malloc_works" != "Xyes"; then
 #  AC_MSG_ERROR([malloc(...) not working, required by AztecOO])
 #fi
fi

if test "X$ac_cv_header_stdexcept" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find stdexcept.  TSFCoreUtils requires at least one or the other.]])
fi  

if test "X$ac_cv_header_memory" != "Xyes" && test "X$ac_cv_header_memory_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find memory or memory.h.  TSFCoreUtils requires at least one or the other.]])
fi  

if test "X$ac_cv_header_map" != "Xyes" && test "X$ac_cv_header_map_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find map or map.h.  TSFCoreUtils requires at least one or the other.]])
fi  

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find vector or vector.h.  TSFCoreUtils requires at least one or the other.]])
fi  

if test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_strstream" != "Xyes" && test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_strstream" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find sstream, strstream, sstream.h, or strstream.h.  TSFCoreUtils requires at least one of these.]])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
AC_CXX_NEW_FOR_SCOPING
#AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

### Include Blas/Lapack example here. ###  For now, developers should see
### The configure.ac file for epetra (Trilinos/packages/epetra) if blas and
### or Lapack are needed.

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
  case $target in
    alpha*)
      if test $CXX = cxx; then
        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
      fi
    ;;
  esac
fi

AC_MSG_CHECKING([for template-specific object files to be linked])
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 AC_MSG_RESULT(none)
else
 AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
fi

AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
##
#  #np# - Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile 
		src/Makefile
		test/Makefile
		test/unary_rtop/Makefile
		test/RTOpT/Makefile
		])

AC_OUTPUT()

# Bye World!
#  #np# change package name
echo "----------------------------------------"
echo "Finished Running TSFCoreUtils Configure Script"
echo "----------------------------------------"
