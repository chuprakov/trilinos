#@HEADER
#@HEADER

EXEEXT = .exe

if BUILD_TEUCHOS_EXPERIMENTAL
if HAVE_TSFCORE_EPETRA
EPETRA_TEST = epetra_adapters
endif
endif

noinst_PROGRAMS = $(EPETRA_TEST)

epetra_adapters_SOURCES = $(srcdir)/cxx_main.cpp

epetra_adapters_DEPENDENCIES = \
  $(top_builddir)/src/libTSFCore.a

epetra_adapters_LDFLAGS =

epetra_adapters_CXXFLAGS = \
  -I$(top_builddir)/src \
  -I$(top_srcdir)/src \
  -I$(top_srcdir)/src/interfaces/Core \
  -I$(top_srcdir)/src/interfaces/Solvers \
  -I$(top_srcdir)/src/interfaces/Nonlin \
  -I$(top_srcdir)/src/utilities/Core \
  -I$(top_srcdir)/src/utilities/Solvers \
  -I$(top_srcdir)/src/utilities/Nonlin \
  -I$(top_srcdir)/src/adapters/mpi-base \
  -I$(top_srcdir)/src/adapters/epetra \
  -I$(top_srcdir)/example/Core \
  -I$(top_srcdir)/example/Nonlin/NP2DSim \
  -I$(top_srcdir)/example/Nonlin/NP4DOpt \
  -I$(top_builddir)/../TSFCoreUtils/src \
  -I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
  -I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
  -I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
  -I$(top_srcdir)/../TSFCoreUtils/src \
  -I$(top_builddir)/../epetra/src \
  -I$(top_srcdir)/../epetra/src \
  -I$(top_builddir)/../teuchos/src \
  -I$(top_srcdir)/../teuchos/src

epetra_adapters_LDADD = \
  $(top_builddir)/src/libTSFCore.a \
  $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(top_builddir)/../teuchos/src/libteuchos.a \
  $(BLAS_LIBS) $(FLIBS)
