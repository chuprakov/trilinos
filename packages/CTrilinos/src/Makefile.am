# @HEADER
# ************************************************************************
# 
# RTOp: Interfaces and Support Software for Vector Reduction Transformation
#       Operations
#                Copyright (2006) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package ctrilinos subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#

if HAVE_MPI
  MPI_SRC_H = \
    $(srcdir)/FEpetra_MpiComm.h

  MPI_SRC_HPP = \
    $(srcdir)/FEpetra_MpiComm_Cpp.hpp

  MPI_SRC = \
    $(srcdir)/FEpetra_MpiComm.cpp
else
  MPI_SRC_H =

  MPI_SRC_HPP =

  MPI_SRC =
endif

CTRILINOS_SRC_H = \
  $(srcdir)/FEpetra_Map.h \
  $(srcdir)/FEpetra_Vector.h \
  $(MPI_SRC_H) \
  $(srcdir)/FEpetra_SerialComm.h

CTRILINOS_SRC_HPP = \
  $(srcdir)/FEpetra_Map_Cpp.hpp \
  $(srcdir)/FEpetra_Vector_Cpp.hpp \
  $(MPI_SRC_HPP) \
  $(srcdir)/FEpetra_SerialComm_Cpp.hpp

CTRILINOS_SRC = \
  $(srcdir)/FEpetra_Map.cpp \
  $(srcdir)/FEpetra_Vector.cpp \
  $(MPI_SRC) \
  $(srcdir)/FEpetra_SerialComm.cpp


#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.CTrilinos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(CTRILINOS_INCLUDES))
else
EXPORT_INCLUDES = $(CTRILINOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Ctrilinos library specifications
# --------------------------------------------------------------------

CTRILINOS_LIB = libctrilinos.a

CTRILINOS_H = \
  $(srcdir)/CTrilinos_ConfigDefs.h \
  $(CTRILINOS_SRC_HPP) \
  $(CTRILINOS_SRC_H)

libctrilinos_a_SOURCES = \
  $(CTRILINOS_SRC)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libctrilinos_a_AR = $(ALTERNATE_AR)

else

libctrilinos_a_AR = $(AR) $(ARFLAGS)

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libctrilinos_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(CTRILINOS_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(CTRILINOS_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: CtrilinosConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = CTrilinos_Config.h

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 


EXTRA_DIST = CTrilinos_version.h CTrilinos_version.c
