# @HEADER
# ************************************************************************
# 
#               TSFCore: Trilinos Solver Framework Core
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
CORE_H = \
	$(srcdir)/Core/TSFCoreAssertOp.hpp\
	$(srcdir)/Core/TSFCoreAssertOpDecl.hpp\
	$(srcdir)/Core/TSFCoreLinearOp.hpp\
	$(srcdir)/Core/TSFCoreLinearOpDecl.hpp\
	$(srcdir)/Core/TSFCoreMultiVector.hpp\
	$(srcdir)/Core/TSFCoreMultiVectorCols.hpp\
	$(srcdir)/Core/TSFCoreMultiVectorColsDecl.hpp\
	$(srcdir)/Core/TSFCoreMultiVectorDecl.hpp\
	$(srcdir)/Core/TSFCoreMultiVectorStdOps.hpp\
	$(srcdir)/Core/TSFCoreMultiVectorStdOpsDecl.hpp\
	$(srcdir)/Core/TSFCoreOpBase.hpp\
	$(srcdir)/Core/TSFCoreOpBaseDecl.hpp\
	$(srcdir)/Core/TSFCoreSerialVector.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorBase.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorBaseDecl.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorDecl.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorSpace.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorSpaceBase.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorSpaceBaseDecl.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorSpaceDecl.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorSpaceFactory.hpp\
	$(srcdir)/Core/TSFCoreSerialVectorSpaceFactoryDecl.hpp\
	$(srcdir)/Core/TSFCoreTypes.hpp\
	$(srcdir)/Core/TSFCoreVector.hpp\
	$(srcdir)/Core/TSFCoreVectorDecl.hpp\
	$(srcdir)/Core/TSFCoreVectorSpace.hpp\
	$(srcdir)/Core/TSFCoreVectorSpaceDecl.hpp\
	$(srcdir)/Core/TSFCoreVectorSpaceFactory.hpp\
	$(srcdir)/Core/TSFCoreVectorSpaceFactoryDecl.hpp\
	$(srcdir)/Core/TSFCoreVectorStdOps.hpp\
	$(srcdir)/Core/TSFCoreVectorStdOpsDecl.hpp\
	$(srcdir)/Core/TSFCore_apply_op_helper.hpp\
	$(srcdir)/Core/TSFCore_apply_op_helper_decl.hpp


CORE = 

SOLVERS_H = \
	$(srcdir)/Solvers/TSFCoreSolversConvergenceTesterDecl.hpp\
	$(srcdir)/Solvers/TSFCoreSolversConvergenceTester.hpp\
	$(srcdir)/Solvers/TSFCoreSolversIterativeLinearSolverDecl.hpp\
	$(srcdir)/Solvers/TSFCoreSolversIterativeLinearSolver.hpp\
	$(srcdir)/Solvers/TSFCoreSolversNormDecl.hpp\
	$(srcdir)/Solvers/TSFCoreSolversNorm.hpp\
	$(srcdir)/Solvers/TSFCoreSolversSolverStateDecl.hpp\
	$(srcdir)/Solvers/TSFCoreSolversSolverState.hpp\
	$(srcdir)/Solvers/TSFCoreSolversTypes.hpp

NONLIN_H = \
	$(srcdir)/Nonlin/TSFCoreNonlinLinearOpWithSolveDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinLinearOpWithSolve.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinLinearSolveOpDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinLinearSolveOp.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrder.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblem.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinTypes.hpp

#Specify locations of include files within Trilinos, but outside TSFCore
### Make sure all of these paths are still necessary
#  
AM_CPPFLAGS = -I$(top_srcdir)/interfaces/Core \
-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack\
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps\
-I$(top_builddir)/../TSFCoreUtils/src/MoochoUtilities

# ------------------------------------------------------------------------
# TSFCore library specifications
# ------------------------------------------------------------------------

#TSFCORE_LIB = libTSFCore.a


#TSFCORE_H = \
#	$(CORE_H) 

#libTSFCore_a_SOURCES = \
#	$(CORE) 



#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#libTSFCore_a_AR = $(ALTERNATE_AR)
else

#libTSFCore_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#libTSFCore_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#lib_LIBRARIES = $(TSFCORE_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#include_HEADERS = $(TSFCORE_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: TSFCoreConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = TSFCoreConfig.h $(CORE_H) $(SOLVERS_H) \
	$(NONLIN_H)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






