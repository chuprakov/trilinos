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

INTERFACES_CORE_H = \
	$(srcdir)/interfaces/Core/TSFCoreAssertOp.hpp\
	$(srcdir)/interfaces/Core/TSFCoreAssertOpDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreLinearOp.hpp\
	$(srcdir)/interfaces/Core/TSFCoreLinearOpDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreMultiVector.hpp\
	$(srcdir)/interfaces/Core/TSFCoreMultiVectorCols.hpp\
	$(srcdir)/interfaces/Core/TSFCoreMultiVectorColsDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreMultiVectorDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreMultiVectorStdOps.hpp\
	$(srcdir)/interfaces/Core/TSFCoreMultiVectorStdOpsDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreOpBase.hpp\
	$(srcdir)/interfaces/Core/TSFCoreOpBaseDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreDotProd.hpp\
	$(srcdir)/interfaces/Core/TSFCoreDotProdDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreScalarProd.hpp\
	$(srcdir)/interfaces/Core/TSFCoreScalarProdDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVector.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorBase.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorBaseDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorSpace.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorSpaceBase.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorSpaceBaseDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorSpaceDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorSpaceFactory.hpp\
	$(srcdir)/interfaces/Core/TSFCoreSerialVectorSpaceFactoryDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreTypes.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVector.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorMultiVector.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorMultiVectorDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceStdBase.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceStdBaseDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpace.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceFactory.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceFactoryDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorStdOps.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorStdOpsDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCore_apply_op_helper.hpp\
	$(srcdir)/interfaces/Core/TSFCore_apply_op_helper_decl.hpp

UTILITIES_CORE_H = \
	$(srcdir)/utilities/Core/TSFCoreExplicitMultiVectorView.hpp\
	$(srcdir)/utilities/Core/TSFCoreExplicitVectorView.hpp\
	$(srcdir)/utilities/Core/TSFCoreMultiVectorAllocator.hpp\
	$(srcdir)/utilities/Core/TSFCoreLinearOpTesterDecl.hpp\
	$(srcdir)/utilities/Core/TSFCoreLinearOpTester.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductVectorSpaceBase.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductMultiVectorBase.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductVectorBase.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductVectorSpace.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductVectorSpaceDecl.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductVector.hpp\
	$(srcdir)/utilities/Core/TSFCoreProductVectorDecl.hpp\
	$(srcdir)/utilities/Core/TSFCoreTestingToolsDecl.hpp\
	$(srcdir)/utilities/Core/TSFCoreTestingTools.hpp

CORE_H = $(INTERFACES_CORE_H) $(UTILITIES_CORE_H)

CORE =

INTERFACES_SOLVERS_H = \
	$(srcdir)/interfaces/Solvers/TSFCoreSolversConvergenceTesterDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversConvergenceTester.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversIterativeLinearSolverDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversIterativeLinearSolver.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversSolverStateDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversSolverState.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversTypes.hpp

UTILITIES_SOLVERS_H = \
	$(srcdir)/utilities/Solvers/TSFCoreSolversAttachConvergenceTesterBaseDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversAttachConvergenceTesterBase.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversBiCGSolverDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversBiCGSolver.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversCGSolverDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversCGSolver.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversGMRESSolverDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversGMRESSolver.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversSimpleGMRESSolverDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversSimpleGMRESSolver.hpp\
	$(srcdir)/../example/Solvers/TSFCoreSolversGmresSolver.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversNormedConvergenceTesterDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversNormedConvergenceTester.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversSummaryOutputterDecl.hpp\
	$(srcdir)/utilities/Solvers/TSFCoreSolversSummaryOutputter.hpp

SOLVERS_H = $(INTERFACES_SOLVERS_H) $(UTILITIES_SOLVERS_H)

SOLVERS =

INTERFACES_NONLIN_H = \
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinLinearOpWithSolveDecl.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinLinearOpWithSolve.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinLinearSolveOpDecl.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinLinearSolveOp.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinNonlinearProblemDecl.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderDecl.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrder.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinNonlinearProblem.hpp\
	$(srcdir)/interfaces/Nonlin/TSFCoreNonlinTypes.hpp

UTILITIES_NONLIN_H = \
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinLinearOpWithSolveIterDecl.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinLinearOpWithSolveIter.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderTesterDecl.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderTester.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinNonlinearProblemTesterDecl.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinNonlinearProblemTester.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinSimpleNewtonSolverDecl.hpp\
	$(srcdir)/utilities/Nonlin/TSFCoreNonlinSimpleNewtonSolver.hpp

NONLIN_H = $(INTERFACES_NONLIN_H) $(UTILITIES_NONLIN_H)

NONLIN =

EXTRA_ADAPTERS_MPI_H = \
	$(srcdir)/adapters/mpi-base/TSFCoreMPIMultiVectorBaseDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIMultiVectorBase.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIMultiVectorStdDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIMultiVectorStd.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorBaseDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorBase.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceBaseDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceBase.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceStdDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceStd.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceFactoryStdDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceFactoryStd.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMultiVectorSerializationDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMultiVectorSerialization.hpp

EXTRA_ADAPTERS_EPETRA_H = \
	$(srcdir)/adapters/epetra/TSFCoreEpetraTypes.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraLinearOp.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraMultiVector.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVector.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVectorSpace.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVectorSpaceFactory.hpp\
	$(srcdir)/adapters/epetra/TSFCore_get_Epetra_MultiVector.hpp

EXTRA_ADAPTERS_EPETRA = \
	$(srcdir)/adapters/epetra/TSFCoreEpetraLinearOp.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraMultiVector.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVector.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVectorSpace.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVectorSpaceFactory.cpp\
	$(srcdir)/adapters/epetra/TSFCore_get_Epetra_MultiVector.cpp

EXTRA_ADAPTERS_AZTECOO_H = \
	$(srcdir)/adapters/aztecoo/Epetra_DiagonalOperator.hpp \
	$(srcdir)/adapters/aztecoo/Epetra_ProductOperator.hpp \
	$(srcdir)/adapters/aztecoo/Epetra_LinearSystemScaler.hpp \
	$(srcdir)/adapters/aztecoo/Epetra_NonlinearProblem.hpp \
	$(srcdir)/adapters/aztecoo/Epetra_NonlinearProblemFirstOrder.hpp \
	$(srcdir)/adapters/aztecoo/Ifpack_PrecGenerator.hpp \
	$(srcdir)/adapters/aztecoo/TSFCoreNonlinEpetraNPFO.hpp \
	$(srcdir)/adapters/aztecoo/TSFCoreNonlinLinearOpWithSolveAztecOO.hpp

EXTRA_ADAPTERS_AZTECOO = \
	$(srcdir)/adapters/aztecoo/Epetra_DiagonalOperator.cpp \
	$(srcdir)/adapters/aztecoo/Epetra_ProductOperator.cpp \
	$(srcdir)/adapters/aztecoo/Epetra_LinearSystemScaler.cpp \
	$(srcdir)/adapters/aztecoo/Epetra_NonlinearProblem.cpp \
	$(srcdir)/adapters/aztecoo/Epetra_NonlinearProblemFirstOrder.cpp \
	$(srcdir)/adapters/aztecoo/Ifpack_PrecGenerator.cpp \
	$(srcdir)/adapters/aztecoo/TSFCoreNonlinEpetraNPFO.cpp \
	$(srcdir)/adapters/aztecoo/TSFCoreNonlinLinearOpWithSolveAztecOO.cpp

if HAVE_TSFCORE_MPI

ADAPTERS_MPI_H = $(EXTRA_ADAPTERS_MPI_H)

else

ADAPTERS_MPI_H = 

endif

if HAVE_TSFCORE_EPETRA

ADAPTERS_EPETRA_H = $(EXTRA_ADAPTERS_EPETRA_H)
ADAPTERS_EPETRA = $(EXTRA_ADAPTERS_EPETRA)

else

ADAPTERS_EPETRA_H =
ADAPTERS_EPETRA =

endif

if HAVE_TSFCORE_AZTECOO

ADAPTERS_AZTECOO_H = $(EXTRA_ADAPTERS_AZTECOO_H)
ADAPTERS_AZTECOO = $(EXTRA_ADAPTERS_AZTECOO)

else

ADAPTERS_AZTECOO_H =
ADAPTERS_AZTECOO =

endif

#Specify locations of include files within Trilinos, but outside TSFCore

AM_CPPFLAGS = \
-I$(srcdir) \
-I$(srcdir)/interfaces/Core \
-I$(srcdir)/interfaces/Solvers \
-I$(srcdir)/interfaces/Nonlin \
-I$(srcdir)/utilities/Core \
-I$(srcdir)/adapters/mpi-base \
-I$(srcdir)/adapters/epetra \
-I$(srcdir)/adapters/aztecoo \
-I$(top_srcdir)/../TSFCoreUtils/src \
-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
-I$(top_builddir)/../TSFCoreUtils/src\
-I$(top_srcdir)/../epetra/src \
-I$(top_builddir)/../epetra/src \
-I$(top_srcdir)/../ifpack/src \
-I$(top_builddir)/../ifpack/src \
-I$(top_srcdir)/../aztecoo/src \
-I$(top_builddir)/../aztecoo/src \
-I$(top_srcdir)/../teuchos \
-I$(top_srcdir)/../teuchos/src \
-I$(top_builddir)/../teuchos/src

# Note, epetra and aztecoo includes only needed and used
# if support for epetra or aztecoo is specifified

# --------------------------------------------------------------------
# TSFCore library specifications
# --------------------------------------------------------------------

TSFCORE_LIB = libTSFCore.a

TSFCORE_H = \
	$(srcdir)/TSFCore_ConfigDefs.hpp \
	$(CORE_H) \
	$(SOLVERS_H) \
	$(NONLIN_H) \
	$(ADAPTERS_MPI_H) \
	$(ADAPTERS_EPETRA_H)\
  $(ADAPTERS_AZTECOO_H)

libTSFCore_a_SOURCES = \
	$(ADAPTERS_EPETRA)\
  $(ADAPTERS_AZTECOO)

EXTRA_libTSFCore_a_SOURCES = \
	$(EXTRA_ADAPTERS_MPI_H) \
	$(EXTRA_ADAPTERS_EPETRA_H) \
	$(EXTRA_ADAPTERS_EPETRA) \
	$(EXTRA_ADAPTERS_AZTECOO_H) \
	$(EXTRA_ADAPTERS_AZTECOO)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libTSFCore_a_AR = $(ALTERNATE_AR)
else

libTSFCore_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

if HAVE_TSFCORE_TTO_ADAPTER

 libTSFCore_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

else 

 libTSFCore_a_LIBADD = $(XTRALDADD)

endif
# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(TSFCORE_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(TSFCORE_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: TSFCoreConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = TSFCoreConfig.h $(TSFCORE_H)

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 


