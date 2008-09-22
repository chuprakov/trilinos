# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                      Didasko Tutorial Package
#                 Copyright (2005) Sandia Corporation
# 
#
# Copyright (2001,2002,2003) Sandia Corportation. Under the terms of Contract 
# DE-AC04-94AL85000, there is a non-exclusive license for use of this 
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
# 
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# 
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# package name, version number, and e-mail address below
# From Trilinos 7.0, the version number should be as for Trilinos itself,
# since the tutorial is so stricly related to Trilinos. Remember to change it 
# a few lines below as well.
AC_INIT(didasko, 9.0, trilinos-help@software.sandia.gov)

# Hello World!
#
echo "----------------------------------------"
echo "Running Didasko Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/Didasko_ConfigDefs.h])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(src/Didasko_config.h:src/Didasko_config.h.in)

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

# Package name and version below
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

#  can eliminate compiler checks below if your package does not use the
#  language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
#AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for extra flags
# You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS(libs, LIBS)
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# Customize the checks for headers needed in package being added.
# NOTE: This section still needs to be updated for didasko, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h. Didasko requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h. Didasko requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h. Didasko requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h. Didasko requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h. Didasko requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h. Didasko requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h. Didasko requires at least one or the other.]])
  fi
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

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
# #np# newp_swahili is an example of a package feature which requires both
# a cpp #define command and one or more source or header files  - ie 
# conditional compilation which requires a change to the makefile.
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in New_Package_config.h and add source 
# and header files as specified in Makefile.am 
#
# To add an option to your package, you will need to modify this file
# (i.e. configure.ac), and src/Makefile.am.  A bootstrap is also required.
# Look for SWAHILI in both files listed above.
#
# TAC_ARG_ENABLE_OPTION performs several tasks, including:
#   1)  Setting a help string
#   2)  Printing a message to standard out
#   3)  Defining HAVE_NEWP_SWAHILI (For use in New_Package_config.h)
#
# AM_CONDITIONAL makes HAVE_NEWP_SWAHILI available for use in Makefile.am
# ------------------------------------------------------------------------
#

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#
#  #np# - These options are specific to Epetra and are left only as examples.
#         If uncommented, the following lines would define HAVE_FATAL_MESSAGES
#         in the (package)_config.h file.  Developers should use a 
#         more specific name such as 'new_package-fatal_messages', rather
#         than 'fatal_messages'.
#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)
#  #np# - These are specific to Epetra and are left only as examples.
#         If uncommented, the following line would define HAVE_ZOLTAN in
#         the (package)_config.h file.  Developers should use a more
#         specific name such as 'ml-zoltan', rather than 'zoltan'.
# TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
# TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)
# TAC_ARG_WITH_PACKAGE(esi, [Enable ESI interface support], ESI, yes)
# TAC_ARG_WITH_PACKAGE(superlu, [Enable SuperLU interface support], SUPERLU, no)

# ------------------------------------------------------------------------
# Interfaces to other Trilinos packages
# ------------------------------------------------------------------------
#
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, epetra, DIDASKO_EPETRA,
  HAVE_DIDASKO_EPETRA, yes,
  [Enables EPETRA examples.], 
  [Enables EPETRA examples.  Can be overridden with --disable-didasko-epetra.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, epetraext, DIDASKO_EPETRAEXT,
  HAVE_DIDASKO_EPETRAEXT, yes,
  [Enables EPETRAEXT examples.],
  [Enables EPETRAEXT examples.  Can be overridden with --disable-didasko-epetraext.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, belos, DIDASKO_BELOS,
  HAVE_DIDASKO_BELOS, yes,
  [Enables BELOS examples.],
  [Enables BELOS examples.  Can be overridden with --disable-didasko-belos.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, nox, DIDASKO_NOX,
  HAVE_DIDASKO_NOX, yes,
  [Enables NOX examples, requires --enable-epetra --enable-nox-epetra],
  [Enables NOX examples, requires --enable-epetra --enable-nox-epetra.  Can be overridden with --disable-didasko-nox.])
#
TAC_ARG_ENABLE_FEATURE_SUB(nox, epetra, [Required when enabling NOX examples], 
  NOX_EPETRA, no)
AM_CONDITIONAL(HAVE_NOX_EPETRA, test "X$ac_cv_use_nox_epetra" = "Xyes")
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, ml, DIDASKO_ML,
  HAVE_DIDASKO_ML, yes,
  [Enables ML examples.],
  [Enables ML examples.  Can be overridden with --disable-didasko-ml.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, amesos, DIDASKO_AMESOS,
  HAVE_DIDASKO_AMESOS, yes,
  [Enables AMESOS examples.],
  [Enables AMESOS examples.  Can be overridden with --disable-didasko-amesos.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, ifpack, DIDASKO_IFPACK,
  HAVE_DIDASKO_IFPACK, yes,
  [Enables IFPACK examples.],
  [Enables IFPACK examples.  Can be overridden with --disable-didasko-ifpack.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, teuchos, DIDASKO_TEUCHOS,
  HAVE_DIDASKO_TEUCHOS, yes,
  [Enables TEUCHOS examples, requires --enable_epetra.],
  [Enables TEUCHOS examples, requires --enable_epetra.  Can be overridden with --disable-didasko-teuchos.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, anasazi, DIDASKO_ANASAZI,
  HAVE_DIDASKO_ANASAZI, yes,
  [Enables ANASAZI examples.],
  [Enables ANASAZI examples.  Can be overridden with --disable-didasko-anasazi.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, belos, DIDASKO_BELOS,
  HAVE_DIDASKO_BELOS, yes,
  [Enables BELOS examples.],
  [Enables BELOS examples.  Can be overridden with --disable-didasko-belos.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, triutils, DIDASKO_TRIUTILS,
  HAVE_DIDASKO_TRIUTILS, yes,
  [Enables TRIUTILS examples.],
  [Enables TRIUTILS examples.  Can be overridden with --disable-didasko-triutils.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, aztecoo, DIDASKO_AZTECOO,
  HAVE_DIDASKO_AZTECOO, yes,
  [Enables AZTECOO examples, requires --enable-teuchos.],
  [Enables AZTECOO examples, requires --enable-teuchos.  Can be overridden with --disable-didasko-aztecoo.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(didasko, tpetra, DIDASKO_TPETRA,
  HAVE_DIDASKO_TPETRA, yes,
  [Enables TPETRA examples.],
  [Enables TPETRA examples.  Can be overridden with --disable-didasko-tpetra.])
#
TAC_ARG_ENABLE_OPTION(latex, [Enables compilation of LaTeX source files.],
		      DIDASKO_LATEX, no)
AM_CONDITIONAL(HAVE_DIDASKO_LATEX, test "X$ac_cv_use_latex" = "Xyes")

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, zoltan, [Enable EpetraExt interface
						support for Zoltan], ZOLTAN,
						no)
AM_CONDITIONAL(BUILD_ZOLTAN, [test "X$ac_cv_use_epetraext_zoltan" != "Xno"])

#
# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( didasko, examples, [Make didasko examples buildable with 'make examples'],
                                 DIDASKO_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_didasko_examples" != "Xno")

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( didasko, tests, [Make didasko tests buildable with 'make tests'],
                                 DIDASKO_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_didasko_tests" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines.
#AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_tests" != "Xno" || test "X$ac_cv_use_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  Only uncomment the following lines if the package
#  creates a subdirectory called cxx_repository for
#  template object files on DEC.  (HINT: try to build on 
#  DEC/CPLANT and see if template-related errors occur)

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
#  Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile 
		Makefile.export
		Makefile.export.didasko 
		examples/Makefile
		examples/amesos/Makefile
		examples/anasazi/Makefile
		examples/aztecoo/Makefile
		examples/belos/Makefile
		examples/epetra/Makefile
		examples/epetraext/Makefile
		examples/ifpack/Makefile
		examples/ml/Makefile
		examples/nox/Makefile
		examples/teuchos/Makefile
		examples/triutils/Makefile
		examples/HBMatrices/Makefile
		examples/tpetra/Makefile
		test/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running Didasko Configure Script"
echo "----------------------------------------"
