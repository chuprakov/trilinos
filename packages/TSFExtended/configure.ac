# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#           TSFExtended: Trilinos Solver Framework Extended
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
AC_INIT(TSFExtended, 1.0, krlong@sandia.gov)

# Hello World!
#  #np# - package name below
echo "----------------------------------------"
echo "Running TSFExtended Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR(/src/Vectors/TSFVectorSpace.hpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/new_package that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #np# - Change file names in next line 
# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

AM_CONFIG_HEADER(src/TSF_config.h:src/TSF_config.h.in)

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

#  #np# - Required version of automake.
AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)



# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the 
#  #np# - option(s) that does (do) not apply.
TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( tsfextended, tests, [Build TSFExtended tests], TSFEXTENDED_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tsfextended_tests" != "Xno")




# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check.  Please note that if you use
#	  F77_FUNC to determine Fortran name mangling, you should not remove
#	  the Fortran compiler check or the check for Fortran flags.  Doing
#	  so will prevent the detection of the proper name mangling in some
#	  cases.
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB



# ------------------------------------------------------------------------
# Checks if ML is enabled
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(TSFExtended, ml, ML, BUILD_ML, yes, [Enable the ML multilevel solvers and preconditioners in TSF], [Enable ML. Can be overridden with --disable-TSFExtended-ml.])




# ------------------------------------------------------------------------
# Get the absolute path to the build directory.  
# 
# We need this to be able to figure out paths to the libraries for the
# other Trilinos packages.
# ------------------------------------------------------------------------

AC_MSG_CHECKING([finding build directory])
if test "X" = "X"; then
  TSF_BUILD_DIRECTORY=`pwd`
  AC_MSG_RESULT($TSF_BUILD_DIRECTORY)
fi


# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)
AC_SUBST(HAVE_ML)
AC_SUBST(TSF_BUILD_DIRECTORY)


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
		Makefile.export.TSFExtended 
		src/Makefile
		test/Makefile
		test/VectorTests/Makefile
		test/Double/Makefile
		test/RectMatrixTest/Makefile
		test/PoissonTest/Makefile
		test/Nonlinear/Makefile
		test/BlockTests/Makefile
  doc/Makefile  
  doc/Doxygen/Makefile
  doc/Doxygen/Doxyfile.user
  doc/Doxygen/Doxyfile.developer
		])

AC_OUTPUT()

# Bye World!
#  #np# change package name
echo "----------------------------------------"
echo "Finished TSFExtended Configure Script"
echo "----------------------------------------"
