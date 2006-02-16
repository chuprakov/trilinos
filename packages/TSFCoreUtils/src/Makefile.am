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

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package tsfcoreutils subdirectory src built successfully."
	@echo ""

MoochoUtilities_SOURCES_H=\
$(srcdir)/MoochoUtilities/Range1D.hpp\
$(srcdir)/MoochoUtilities/StandardAggregationMacros.hpp\
$(srcdir)/MoochoUtilities/StandardCompositionRelationshipsPack.hpp\
$(srcdir)/MoochoUtilities/ReleaseResource_ref_count_ptr.hpp\
$(srcdir)/MoochoUtilities/ReleaseResource.hpp \
$(srcdir)/MoochoUtilities/ObjectDB.hpp

MoochoUtilities_SOURCES=\
$(srcdir)/MoochoUtilities/Range1D.cpp\
$(srcdir)/MoochoUtilities/ObjectDB.cpp

RTOp_C_SOURCES_H=\
$(srcdir)/RTOpPack/RTOp_config.h\
$(srcdir)/RTOpPack/RTOp_parallel_helpers.h\
$(srcdir)/RTOpPack/RTOp_mpi.h

RTOp_C_SOURCES=\
$(srcdir)/RTOpPack/RTOp_parallel_helpers.c\
$(srcdir)/RTOpPack/RTOp_mpi.c

RTOp_CXX_SOURCES_H=\
$(srcdir)/RTOpPack/RTOpFactoryTmpl.hpp\
$(srcdir)/RTOpPack/RTOpPackTypes.hpp\
$(srcdir)/RTOpPack/print_sub_vector.hpp\
$(srcdir)/RTOpPack/RTOpServerDecl.hpp\
$(srcdir)/RTOpPack/RTOpServer.hpp\
$(srcdir)/RTOpPack/RTOpUnaryFuncPtrDecl.hpp\
$(srcdir)/RTOpPack/RTOpUnaryFuncPtr.hpp\
$(srcdir)/RTOpPack/RTOpPack_RTOpTDecl.hpp\
$(srcdir)/RTOpPack/RTOpPack_RTOpT.hpp\
$(srcdir)/RTOpPack/RTOpPack_RTOpTHelpers.hpp\
$(srcdir)/RTOpPack/RTOpPack_MPI_apply_op_decl.hpp\
$(srcdir)/RTOpPack/RTOpPack_MPI_apply_op.hpp

RTOp_CXX_SOURCES=\
$(srcdir)/RTOpPack/print_sub_vector.cpp\
$(srcdir)/RTOpPack/RTOpPack_MPI_apply_op.cpp

RTOp_STDOPS_SOURCES_CPP=

RTOp_STDOPS_SOURCES_HPP=\
$(srcdir)/RTOpStdOps/RTOpPack_ROpDotProd.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpGetSubVector.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpMax.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpMaxIndex.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpMaxIndexLessThanBound.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpMin.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpMinIndex.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpMinIndexGreaterThanBound.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpNorm1.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpNorm2.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpNormInf.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpSum.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_ROpWeightedNorm2.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpAXPY.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpAbs.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpAddScalar.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpAssignVectors.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpEleWiseDivide.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpEleWiseProd.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpLinearCombination.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpRandomize.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpReciprocal.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpAssignScalar.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpScaleVector.hpp\
$(srcdir)/RTOpStdOps/RTOpPack_TOpSetSubVector.hpp

AM_CPPFLAGS = -I$(top_builddir)/src \
  -I$(top_srcdir)/src/MoochoUtilities \
  -I$(top_srcdir)/src/RTOpPack \
  -I$(top_srcdir)/src/RTOpStdOps \
  -I$(top_builddir)/../teuchos/src \
  -I$(top_srcdir)/../teuchos \
  -I$(top_srcdir)/../teuchos/src

SOURCES_H=\
  TSFCoreUtils_ConfigDefs.hpp\
  $(MoochoUtilities_SOURCES_H)\
  $(RTOp_C_SOURCES_H)\
  $(RTOp_CXX_SOURCES_H)\
	$(RTOp_STDOPS_SOURCES_H)\
	$(RTOp_STDOPS_SOURCES_HPP)

TSFCOREUTILS_LIB = libTSFCoreUtils.a

libTSFCoreUtils_a_SOURCES = \
  $(MoochoUtilities_SOURCES)\
  $(RTOp_C_SOURCES)\
  $(RTOp_CXX_SOURCES)\
  $(RTOp_STDOPS_SOURCES)\
  $(RTOp_STDOPS_SOURCES_CPP)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libTSFCoreUtils_a_AR = $(ALTERNATE_AR)
else

libTSFCoreUtils_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libTSFCoreUtils_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(TSFCOREUTILS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(SOURCES_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: TSFCoreUtils_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = TSFCoreUtilsConfig.h $(SOURCES_H)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

