#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = rcp_test

rcp_test_SOURCES = $(srcdir)/cxx_main.cpp

rcp_test_DEPENDENCIES = $(top_builddir)/src/libTSFCoreUtils.a

rcp_test_CXXFLAGS = -I$(top_srcdir)/src/MoochoUtilities

rcp_test_LDADD = -L$(top_builddir)/src -lTSFCoreUtils $(FLIBS)


