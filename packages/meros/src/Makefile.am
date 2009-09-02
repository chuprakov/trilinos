# @HEADER
# ************************************************************************
# 
#              Meros: Segregated Preconditioning Package
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

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package meros subdirectory src-thrya built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside meros. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#
CORE_H = \
	$(srcdir)/Meros_AztecSolveStrategy.hpp        \
	$(srcdir)/Meros_ConfigDefs.h                  \
	$(srcdir)/Meros_IdentityOperator.hpp          \
	$(srcdir)/Meros_InverseOperator.hpp           \
	$(srcdir)/Meros_LinearSolver.hpp              \
	$(srcdir)/Meros_LinearSolveStrategy.hpp       \
	$(srcdir)/Meros_LSCOperatorSource.h           \
	$(srcdir)/Meros_LSCPreconditionerFactory.h    \
	$(srcdir)/Meros_PCDOperatorSource.h           \
	$(srcdir)/Meros_PCDPreconditionerFactory.h    \
	$(srcdir)/Meros_Preconditioner.hpp	      \
	$(srcdir)/Meros_PreconditioningStrategy.hpp   \
	$(srcdir)/Meros_SIMPLEOperatorSource.h        \
	$(srcdir)/Meros_SIMPLEPreconditionerFactory.h \
	$(srcdir)/Meros_ZeroOperator.hpp              

CORE = \
	$(srcdir)/Meros_LSCOperatorSource.cpp         \
	$(srcdir)/Meros_LSCPreconditionerFactory.cpp  \
	$(srcdir)/Meros_PCDOperatorSource.cpp         \
	$(srcdir)/Meros_PCDPreconditionerFactory.cpp  \
	$(srcdir)/Meros_SIMPLEOperatorSource.cpp      \
	$(srcdir)/Meros_SIMPLEPreconditionerFactory.cpp  

UTIL_H = \
	$(srcdir)/Meros_Version.h

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#  Files which are included in select builds are typically have
#  names that start with EXTRA
#
#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif


USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 

include $(top_builddir)/Makefile.export.meros
AM_CPPFLAGS = $(MEROS_INCLUDES)

#AM_CPPFLAGS = -I$(top_builddir)/include
#AM_CPPFLAGS = -I$(srcdir) \
#              -I$(top_builddir)/../meros/src                     \
#              -I$(top_srcdir)/../teuchos/src		         \
#              -I$(top_builddir)/../teuchos/src                   \
#              -I$(top_srcdir)/../rtop/src	                 \
#              -I$(top_builddir)/../rtop/src                      \
#              -I$(top_srcdir)/../rtop/src/interfaces	         \
#              -I$(top_builddir)/../rtop/src/interfaces           \
#              -I$(top_srcdir)/../rtop/src/ops_lib	         \
#              -I$(top_builddir)/../rtop/src/ops_lib              \
#              -I$(top_srcdir)/../rtop/src/support	         \
#              -I$(top_builddir)/../rtop/src/support              \
#              -I$(top_srcdir)/../thyra/src                       \
#              -I$(top_builddir)/../thyra/src/                    \
# -I$(top_builddir)/../thyra/src/interfaces/operator_solve/ana/extended     \
# -I$(top_srcdir)/../thyra/src/interfaces/operator_solve/ana/extended       \
# -I$(top_builddir)/../thyra/src/interfaces/operator_solve/ana/fundamental  \
# -I$(top_srcdir)/../thyra/src/interfaces/operator_solve/ana/fundamental    \
# -I$(top_builddir)/../thyra/src/interfaces/operator_vector/ana/extended    \
# -I$(top_srcdir)/../thyra/src/interfaces/operator_vector/ana/extended      \
# -I$(top_builddir)/../thyra/src/interfaces/operator_vector/ana/fundamental \
# -I$(top_srcdir)/../thyra/src/interfaces/operator_vector/ana/fundamental   \
# -I$(top_builddir)/../thyra/src/interfaces/operator_vector/app             \
# -I$(top_srcdir)/../thyra/src/interfaces/operator_vector/app               \
# -I$(top_builddir)/../thyra/src/support/operator_solve/client_support      \
# -I$(top_srcdir)/../thyra/src/support/operator_solve/client_support        \
# -I$(top_builddir)/../thyra/src/support/operator_vector/adapter_support    \
# -I$(top_srcdir)/../thyra/src/support/operator_vector/adapter_support      \
# -I$(top_builddir)/../thyra/src/support/operator_vector/client_support     \
# -I$(top_srcdir)/../thyra/src/support/operator_vector/client_support       \
#              -I$(top_srcdir)/../epetra/src		  \
#              -I$(top_builddir)/../epetra/src             \
#              -I$(top_srcdir)/../epetra/thyra/src	  \
#              -I$(top_builddir)/../epetra/thyra/src       \
#              -I$(top_srcdir)/../epetraext/src		  \
#              -I$(top_builddir)/../epetraext/src          \
#              -I$(top_srcdir)/../epetraext/src/inout	  \
#              -I$(top_builddir)/../epetraext/inout        \
#              -I$(top_srcdir)/../aztecoo/src              \
#              -I$(top_srcdir)/../ml/src/                  \
#              -I$(top_srcdir)/../ml/src/Coarsen           \
#              -I$(top_srcdir)/../ml/src/Comm              \
#              -I$(top_srcdir)/../ml/src/FEGrid            \
#              -I$(top_srcdir)/../ml/src/Include           \
#              -I$(top_srcdir)/../ml/src/Krylov            \
#              -I$(top_srcdir)/../ml/src/Main              \
#              -I$(top_srcdir)/../ml/src/MatrixFree        \
#              -I$(top_srcdir)/../ml/src/MLAPI             \
#              -I$(top_srcdir)/../ml/src/Mrtr              \
#              -I$(top_srcdir)/../ml/src/NonlinML          \
#              -I$(top_srcdir)/../ml/src/Operator          \
#              -I$(top_srcdir)/../ml/src/Runtime           \
#              -I$(top_srcdir)/../ml/src/Smoother          \
#              -I$(top_srcdir)/../ml/src/Utils             \
#              -I$(top_builddir)/../aztecoo/src              \
#              -I$(top_builddir)/../ml/src/                  \
#              -I$(top_builddir)/../ml/src/Coarsen           \
#              -I$(top_builddir)/../ml/src/Comm              \
#              -I$(top_builddir)/../ml/src/FEGrid            \
#              -I$(top_builddir)/../ml/src/Include           \
#              -I$(top_builddir)/../ml/src/Krylov            \
#              -I$(top_builddir)/../ml/src/Main              \
#              -I$(top_builddir)/../ml/src/MatrixFree        \
#              -I$(top_builddir)/../ml/src/MLAPI             \
#              -I$(top_builddir)/../ml/src/Mrtr              \
#              -I$(top_builddir)/../ml/src/NonlinML          \
#              -I$(top_builddir)/../ml/src/Operator          \
#              -I$(top_builddir)/../ml/src/Runtime           \
#              -I$(top_builddir)/../ml/src/Smoother          \
#              -I$(top_builddir)/../ml/src/Utils             \
#              -I$(top_builddir)/../aztecoo/src              \
#              -I$(top_builddir)/../aztecoo/thyra/src  



# ------------------------------------------------------------------------
# meros  library specifications
# ------------------------------------------------------------------------

MEROS_LIB = libmeros.a


MEROS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H)

libmeros_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI)


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libmeros_a_AR = $(ALTERNATE_AR)
else

libmeros_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libmeros_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(MEROS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(MEROS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: meros_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(includedir)
nodist_execinclude_HEADERS = Meros_Config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






