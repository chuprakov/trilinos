#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = TSFCoreSolversGmresTest

TSFCoreSolversGmresTest_SOURCES = $(srcdir)/TSFCoreSolversGmresTest.cpp

TSFCoreSolversGmresTest_DEPENDENCIES = \
$(top_builddir)/src/libTSFCore.a \
$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
$(top_builddir)/../teuchos/src/libteuchos.a \
$(top_builddir)/../triutils/src/libtriutils.a \
$(top_builddir)/../epetra/src/libepetra.a

TSFCoreSolversGmresTest_CXXFLAGS = \
-I$(top_srcdir)/src \
-I$(top_builddir)/src \
-I$(top_srcdir)/src/interfaces/Core \
-I$(top_srcdir)/src/adapters/epetra \
-I$(top_srcdir)/src/adapters/mpi-base \
-I$(top_srcdir)/../TSFCoreUtils/src \
-I$(top_builddir)/../TSFCoreUtils/src \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
-I$(top_srcdir)/../teuchos/src \
-I$(top_builddir)/../teuchos/src \
-I$(top_srcdir)/../epetra/src \
-I$(top_builddir)/../epetra/src \
-I$(top_srcdir)/../triutils/src \
-I$(top_builddir)/../triutils/src

TSFCoreSolversGmresTest_LDADD = \
	$(top_builddir)/src/libTSFCore.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


