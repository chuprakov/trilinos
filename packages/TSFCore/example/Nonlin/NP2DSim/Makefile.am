#@HEADER
#@HEADER


EXEEXT = .exe

if BUILD_TEUCHOS
NONLIN_EX = TSFCoreNonlinNP2DSimTest
endif

noinst_PROGRAMS = $(NONLIN_EX)

TSFCoreNonlinNP2DSimTest_SOURCES = $(srcdir)/TSFCoreNonlinNP2DSimTestMain.cpp \
	$(srcdir)/TSFCoreNonlinNP2DSimDecl.hpp \
	$(srcdir)/TSFCoreNonlinNP2DSim.hpp

TSFCoreNonlinNP2DSimTest_DEPENDENCIES =\
$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a

TSFCoreNonlinNP2DSimTest_CXXFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/interfaces/Core \
	-I$(top_srcdir)/src/interfaces/Solvers \
	-I$(top_srcdir)/src/interfaces/Nonlin \
	-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
	-I$(top_srcdir)/src/utilities/Nonlin \
	-I$(top_srcdir)/src/utilities/Core \
	-I$(top_srcdir)/src/utilities/Solvers \
	-I$(top_srcdir)/../TSFCoreUtils/src \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
	-I$(top_builddir)/../TSFCoreUtils/src \
    -I$(top_builddir)/../teuchos/src \
    -I$(top_srcdir)/../teuchos \
    -I$(top_srcdir)/../teuchos/src

TSFCoreNonlinNP2DSimTest_LDADD = \
	$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LIBS) $(FLIBS)


