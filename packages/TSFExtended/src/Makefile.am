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
	$(srcdir)/Vectors/TSFLinearCombination.hpp\
	$(srcdir)/Vectors/TSFAccessibleVector.hpp\
	$(srcdir)/Vectors/TSFProductVectorSpace.hpp\
	$(srcdir)/Vectors/TSFVector.hpp \
	$(srcdir)/Vectors/TSFVectorSpace.hpp \
	$(srcdir)/Vectors/TSFVectorType.hpp \
	$(srcdir)/Vectors/TSFVectorTypeExtensions.hpp

VECTORS = \
	$(srcdir)/Vectors/TSFProductVectorSpace.cpp

UTILS_H = \
	$(srcdir)/Utils/TSFConfigDefs.hpp \
	$(srcdir)/Utils/TSFDescribable.hpp \
	$(srcdir)/Utils/TSFDescribableByTypeID.hpp \
	$(srcdir)/Utils/TSFDouble.hpp \
	$(srcdir)/Utils/TSFHandle.hpp \
	$(srcdir)/Utils/TSFHandleable.hpp \
	$(srcdir)/Utils/TSFObjectWithVerbosity.hpp \
	$(srcdir)/Utils/TSFPrintable.hpp

UTILS = 

CONCRETE_H = \
	$(srcdir)/Concrete/TSFEpetraMatrix.hpp\
	$(srcdir)/Concrete/TSFEpetraVector.hpp\
	$(srcdir)/Concrete/TSFEpetraVectorSpace.hpp\
	$(srcdir)/Concrete/TSFEpetraVectorType.hpp \
	$(srcdir)/Concrete/TSFIfpackOperator.hpp

CONCRETE = \
	$(srcdir)/Concrete/TSFEpetraMatrix.cpp\
	$(srcdir)/Concrete/TSFEpetraVector.cpp\
	$(srcdir)/Concrete/TSFEpetraVectorSpace.cpp\
	$(srcdir)/Concrete/TSFEpetraVectorType.cpp \
	$(srcdir)/Concrete/TSFIfpackOperator.cpp

OPERATORS_H = \
	$(srcdir)/Operators/TSFBlockOperator.hpp \
	$(srcdir)/Operators/TSFComposedOperator.hpp \
	$(srcdir)/Operators/TSFDiagonalOperator.hpp\
	$(srcdir)/Operators/TSFExplicitlyTransposeableOp.hpp\
	$(srcdir)/Operators/TSFIdentityOperator.hpp \
	$(srcdir)/Operators/TSFILUFactorizableOp.hpp \
	$(srcdir)/Operators/TSFInverseOperator.hpp \
	$(srcdir)/Operators/TSFLinearOperator.hpp \
	$(srcdir)/Operators/TSFLinearOperatorDecl.hpp \
	$(srcdir)/Operators/TSFLinearOperatorImpl.hpp \
	$(srcdir)/Operators/TSFLoadableMatrix.hpp\
	$(srcdir)/Operators/TSFNonlinearOperator.hpp\
	$(srcdir)/Operators/TSFNonlinearOperatorBase.hpp\
	$(srcdir)/Operators/TSFOpDescribableByTypeID.hpp \
	$(srcdir)/Operators/TSFRowAccessibleOp.hpp \
	$(srcdir)/Operators/TSFScaledOperator.hpp\
	$(srcdir)/Operators/TSFSumOperator.hpp\
	$(srcdir)/Operators/TSFTransposeOperator.hpp\
	$(srcdir)/Operators/TSFZeroOperator.hpp 

OPERATORS = \
	$(srcdir)/Operators/TSFBlockOperator.cpp 

LINEARSOLVERS_H = \
	$(srcdir)/LinearSolvers/TSFAztecSolver.hpp \
	$(srcdir)/LinearSolvers/TSFBICGSTABSolver.hpp \
	$(srcdir)/LinearSolvers/TSFGMRESSolver.hpp \
	$(srcdir)/LinearSolvers/TSFIterativeSolver.hpp \
	$(srcdir)/LinearSolvers/TSFKrylovSolver.hpp \
	$(srcdir)/LinearSolvers/TSFLinearSolver.hpp \
	$(srcdir)/LinearSolvers/TSFLinearSolverBase.hpp \
	$(srcdir)/LinearSolvers/TSFLinearSolverBaseDecl.hpp \
	$(srcdir)/LinearSolvers/TSFLinearSolverBaseImpl.hpp \
	$(srcdir)/LinearSolvers/TSFLinearSolverBuilder.hpp \
	$(srcdir)/LinearSolvers/TSFSolverState.hpp 

LINEARSOLVERS = \
	$(srcdir)/LinearSolvers/TSFAztecSolver.cpp \
	$(srcdir)/LinearSolvers/TSFLinearSolverBuilder.cpp 


PRECONDITIONERS_H = \
	$(srcdir)/Preconditioners/TSFGenericLeftPreconditioner.hpp \
	$(srcdir)/Preconditioners/TSFGenericRightPreconditioner.hpp \
	$(srcdir)/Preconditioners/TSFILUKPreconditionerFactory.hpp \
	$(srcdir)/Preconditioners/TSFPreconditioner.hpp \
	$(srcdir)/Preconditioners/TSFPreconditionerBase.hpp \
	$(srcdir)/Preconditioners/TSFPreconditionerFactory.hpp \
	$(srcdir)/Preconditioners/TSFPreconditionerFactoryBase.hpp 


PRECONDITIONERS = 


#
# Read macros that are common to building source and tests. We define these
# in Makefile.common to avoid maintaining the lists in multiple files.
#
include $(top_srcdir)/Makefile.common



#np# List locations of include files within Trilinos but outside of this 
#np# package
#np# - Add packages that this package is dependent on
#np#   Remove epetra if this package does not need it
#  
AM_CPPFLAGS = $(COMMON_CXXFLAGS)


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
	$(OPERATORS_H) \
	$(LINEARSOLVERS_H) \
	$(PRECONDITIONERS_H) \
	$(CONCRETE_H) 

#np# replace new_package with the name of the package being autotool'ed here
libTSFExtended_a_SOURCES = \
	$(UTILS) \
	$(VECTORS) \
	$(OPERATORS) \
	$(LINEARSOLVERS) \
	$(PRECONDITIONERS) \
	$(CONCRETE) 


#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libTSFExtended_a_SOURCES = 



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
nodist_execinclude_HEADERS = TSF_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






