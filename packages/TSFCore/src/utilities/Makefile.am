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
CORE_H = \
	$(srcdir)/Core/TSFCoreExplicitVectorView.hpp\
	$(srcdir)/Core/TSFCoreMultiVectorAllocator.hpp\
	$(srcdir)/Core/TSFCoreTestingToolsDecl.hpp\
	$(srcdir)/Core/TSFCoreTestingTools.hpp

CORE = 

SOLVERS_H = \
	$(srcdir)/Solvers/TSFCoreSolversBiCGSolverDecl.hpp\
	$(srcdir)/Solvers/TSFCoreSolversBiCGSolver.hpp\
	$(srcdir)/Solvers/TSFCoreSolversNormedConvergenceTesterDecl.hpp\
	$(srcdir)/Solvers/TSFCoreSolversNormedConvergenceTester.hpp

NONLIN_H = \
	$(srcdir)/Nonlin/TSFCoreNonlinLinearOpWithSolveIterDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinLinearOpWithSolveIter.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderTesterDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemFirstOrderTester.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemTesterDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinNonlinearProblemTester.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinSimpleNewtonSolverDecl.hpp\
	$(srcdir)/Nonlin/TSFCoreNonlinSimpleNewtonSolver.hpp

#Specify locations of include files within Trilinos, but outside TSFCore
### Make sure all of these paths are still necessary
#  
AM_CPPFLAGS = -I$(top_srcdir)/interfaces/Core -I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities\
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
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

#pkginclude_HEADERS = $(TSFCORE_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: TSFCoreConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/TSFCore
nodist_pkgexecinclude_HEADERS = $(CORE_H) $(SOLVERS_H) $(NONLIN_H)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






