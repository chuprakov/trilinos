#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = unary_rtop

unary_rtop_SOURCES = $(srcdir)/cxx_main.cpp

unary_rtop_DEPENDENCIES = \
  $(top_builddir)/src/libTSFCoreUtils.a

unary_rtop_LDFLAGS = \
  -L$(top_builddir)/../teuchos/src \
  -L$(top_builddir)/src

unary_rtop_CXXFLAGS = \
  -I$(top_builddir)/src \
  -I$(top_srcdir)/src/RTOpPack \
  -I$(top_srcdir)/src/MoochoUtilities \
  -I$(top_srcdir)/src \
  -I$(top_builddir)/../teuchos/src \
  -I$(top_srcdir)/../teuchos/src

unary_rtop_LDADD = $(unary_rtop_LDFLAGS) -L$(top_builddir)/src -lTSFCoreUtils -lteuchos $(FLIBS)
