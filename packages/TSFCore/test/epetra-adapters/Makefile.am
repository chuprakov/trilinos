#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = epetra_adapters

epetra_adapters_SOURCES = $(srcdir)/cxx_main.cpp

epetra_adapters_DEPENDENCIES = \
  $(top_builddir)/src/libTSFCore.a

epetra_adapters_LDFLAGS = \
  -L$(top_builddir)/../teuchos/src \
  -L$(top_builddir)/../TSFCoreUtils/src \
  -L$(top_builddir)/../epetra/src \
  -L$(top_builddir)/src

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

epetra_adapters_LDADD = $(epetra_adapters_LDFLAGS) -lTSFCore -lTSFCoreUtils -lepetra -lteuchos -lblas $(FLIBS)
