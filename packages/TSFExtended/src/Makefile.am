# @HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
# 
# ************************************************************************
# @HEADER

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.

VECTORS_H = \
	$(srcdir)/Vectors/TSFLoadableVector.hpp\
	$(srcdir)/Vectors/TSFAccessibleVector.hpp\
	$(srcdir)/Vectors/TSFVector.hpp
	$(srcdir)/Vectors/TSFVectorSpace.hpp

VECTORS = 

UTILS_H = \
	$(srcdir)/Utils/TSFConfigDefs.hpp\
	$(srcdir)/Utils/TSFDescribable.hpp\
	$(srcdir)/Utils/TSFDouble.hpp\
	$(srcdir)/Utils/TSFHandle.hpp\
	$(srcdir)/Utils/TSFHandleable.hpp\
	$(srcdir)/Utils/TSFPrintable.hpp

UTILS = 

CONCRETE_H = \
	$(srcdir)/Concrete/TSFEpetraVector.hpp\
	$(srcdir)/Concrete/TSFEpetraVectorSpace.hpp

CONCRETE = \
	$(srcdir)/Concrete/TSFEpetraVector.cpp\
	$(srcdir)/Concrete/TSFEpetraVectorSpace.cpp

TSFEXTENDED_INCLUDES =  \
-I$(top_builddir)/src \
-I$(top_srcdir)/src/Utils \
-I$(top_srcdir)/src/Vectors \
-I$(top_builddir)/src/Utils \
-I$(top_builddir)/src/Vectors 

TSFCORE_INCLUDES =  \
-I$(top_builddir)/../TSFCore/src \
-I$(top_builddir)/../TSFCore/src/interfaces/Core \
-I$(top_builddir)/../TSFCore/src/utilities/Core \
-I$(top_builddir)/../TSFCore/src/adapters/aztecoo \
-I$(top_builddir)/../TSFCore/src/adapters/epetra \
-I$(top_builddir)/../TSFCore/src/adapters/mpi-base \
-I$(top_srcdir)/../TSFCore/src \
-I$(top_srcdir)/../TSFCore/src/interfaces/Core \
-I$(top_srcdir)/../TSFCore/src/utilities/Core \
-I$(top_srcdir)/../TSFCore/src/adapters/aztecoo \
-I$(top_srcdir)/../TSFCore/src/adapters/epetra \
-I$(top_srcdir)/../TSFCore/src/adapters/mpi-base \
-I$(top_srcdir)/../TSFCoreUtils/src \
-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
-I$(top_builddir)/../TSFCoreUtils/src



#np# List locations of include files within Trilinos but outside of this 
#np# package
#np# - Add packages that this package is dependent on
#np#   Remove epetra if this package does not need it
#  
AM_CPPFLAGS = $(TSFEXTENDED_INCLUDES) $(TSFCORE_INCLUDES)\
 -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src\
 -I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# New_Package library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
TSF_EXTENDED_LIB = libTSFExtended.a

#np# replace new_package with the name of the package being autotool'ed here
TSF_EXTENDED_H = \
	$(UTILS_H) \
	$(VECTORS_H) \
	$(CONCRETE_H) 

#np# replace new_package with the name of the package being autotool'ed here
libTSFExtended_a_SOURCES = \
	$(UTILS) \
	$(VECTORS) \
	$(CONCRETE) 


#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libTSFExtended_a_SOURCES = 


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace new_package with the name of the package being autotool'ed here
libTSFExtended_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libTSFExtended_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libTSFExtended_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(TSF_EXTENDED_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(TSF_EXTENDED_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = TSFExtended_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






