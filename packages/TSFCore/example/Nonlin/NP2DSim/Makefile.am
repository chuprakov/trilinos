#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = TSFCoreNonlinNP2DSimTest

TSFCoreNonlinNP2DSimTest_SOURCES = $(srcdir)/TSFCoreNonlinNP2DSimTestMain.cpp \
	$(srcdir)/TSFCoreNonlinNP2DSimDecl.hpp \
	$(srcdir)/TSFCoreNonlinNP2DSim.hpp

TSFCoreNonlinNP2DSimTest_DEPENDENCIES =\
$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a

TSFCoreNonlinNP2DSimTest_CXXFLAGS = -I$(top_srcdir)/src/interfaces/Core \
	-I$(top_srcdir)/src/interfaces/Solvers \
	-I$(top_srcdir)/src/interfaces/Nonlin \
	-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
	-I$(top_srcdir)/src/extended/Nonlin \
	-I$(top_srcdir)/src/extended/Core \
	-I$(top_srcdir)/src/extended/Solvers \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
	-I$(top_builddir)/../TSFCoreUtils/src/MoochoUtilities

TSFCoreNonlinNP2DSimTest_LDADD = \
	$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
	$(LIBS) $(FLIBS)


