# @HEADER
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
	$(srcdir)/interfaces/Core/TSFCoreVectorSpace.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceFactory.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorSpaceFactoryDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorStdOps.hpp\
	$(srcdir)/interfaces/Core/TSFCoreVectorStdOpsDecl.hpp\
	$(srcdir)/interfaces/Core/TSFCore_apply_op_helper.hpp\
	$(srcdir)/interfaces/Core/TSFCore_apply_op_helper_decl.hpp

EXTENDED_CORE_H = \
	$(srcdir)/extended/Core/TSFCoreExplicitVectorView.hpp\
	$(srcdir)/extended/Core/TSFCoreMultiVectorAllocator.hpp\
	$(srcdir)/extended/Core/TSFCoreTestingToolsDecl.hpp\
	$(srcdir)/extended/Core/TSFCoreTestingTools.hpp

CORE_H = $(INTERFACES_CORE_H) $(EXTENDED_CORE_H)

CORE =

INTERFACES_SOLVERS_H = \
	$(srcdir)/interfaces/Solvers/TSFCoreSolversConvergenceTesterDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversConvergenceTester.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversIterativeLinearSolverDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversIterativeLinearSolver.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversNormDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversNorm.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversSolverStateDecl.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversSolverState.hpp\
	$(srcdir)/interfaces/Solvers/TSFCoreSolversTypes.hpp

EXTENDED_SOLVERS_H = \
	$(srcdir)/extended/Solvers/TSFCoreSolversBiCGSolverDecl.hpp\
	$(srcdir)/extended/Solvers/TSFCoreSolversBiCGSolver.hpp\
	$(srcdir)/extended/Solvers/TSFCoreSolversNormedConvergenceTesterDecl.hpp\
	$(srcdir)/extended/Solvers/TSFCoreSolversNormedConvergenceTester.hpp

SOLVERS_H = $(INTERFACES_SOLVERS_H) $(EXTENDED_SOLVERS_H)

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

EXTENDED_NONLIN_H = \
	$(srcdir)/extended/Nonlin/TSFCoreNonlinLinearOpWithSolveIterDecl.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinLinearOpWithSolveIter.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderTesterDecl.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderTester.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinNonlinearProblemTesterDecl.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinNonlinearProblemTester.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinSimpleNewtonSolverDecl.hpp\
	$(srcdir)/extended/Nonlin/TSFCoreNonlinSimpleNewtonSolver.hpp

NONLIN_H = $(INTERFACES_NONLIN_H) $(EXTENDED_NONLIN_H)

NONLIN =

EXTRA_ADAPTERS_MPI_H = \
	$(srcdir)/adapters/mpi-base/TSFCoreMPIMultiVectorBaseDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIMultiVectorBase.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorBaseDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorBase.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceBaseDecl.hpp\
	$(srcdir)/adapters/mpi-base/TSFCoreMPIVectorSpaceBase.hpp

EXTRA_ADAPTERS_EPETRA_H = \
	$(srcdir)/adapters/epetra/TSFCoreEpetraLinearOp.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraMultiVector.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVector.hpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVectorSpace.hpp

EXTRA_ADAPTERS_EPETRA = \
	$(srcdir)/adapters/epetra/TSFCoreEpetraLinearOp.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraMultiVector.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVector.cpp\
	$(srcdir)/adapters/epetra/TSFCoreEpetraVectorSpace.cpp

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

#Specify locations of include files within Trilinos, but outside TSFCore

AM_CPPFLAGS = -I$(srcdir)/interfaces/Core\
-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities\
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack\
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps\
-I$(top_builddir)/../TSFCoreUtils/src/MoochoUtilities

# --------------------------------------------------------------------
# TSFCore library specifications
# --------------------------------------------------------------------

TSFCORE_LIB = libTSFCore.a


TSFCORE_H = \
	$(CORE_H) \
	$(SOLVERS_H) \
	$(NONLIN_H) \
	$(ADAPTERS_MPI_H) \
	$(ADAPTERS_EPETRA_H) 

libTSFCore_a_SOURCES = \
	$(ADAPTERS_EPETRA) 

EXTRA_libTSFCore_a_SOURCES = \
	$(EXTRA_ADAPTERS_MPI_H)
	$(EXTRA_ADAPTERS_EPETRA_H)
	$(EXTRA_ADAPTERS_EPETRA)

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

libTSFCore_a_LIBADD = $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(TSFCORE_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# --------------------------------------------------------------------

pkginclude_HEADERS = $(TSFCORE_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: TSFCoreConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/TSFCore
nodist_pkgexecinclude_HEADERS = TSFCoreConfig.h $(TSFCORE_H)

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 


