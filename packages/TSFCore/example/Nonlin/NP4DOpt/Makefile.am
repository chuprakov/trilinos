#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = TSFCoreNonlinNP4DOptTest

TSFCoreNonlinNP4DOptTest_SOURCES = $(srcdir)/TSFCoreNonlinNP4DOptTestMain.cpp \
	$(srcdir)/TSFCoreNonlinNP4DOptDecl.hpp \
	$(srcdir)/TSFCoreNonlinNP4DOpt.hpp

TSFCoreNonlinNP4DOptTest_DEPENDENCIES =\
$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a

TSFCoreNonlinNP4DOptTest_CXXFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/example/Nonlin/NP2DSim \
	-I$(top_srcdir)/src/interfaces/Core \
	-I$(top_srcdir)/src/interfaces/Solvers \
	-I$(top_srcdir)/src/interfaces/Nonlin \
	-I$(top_srcdir)/src/utilities/Nonlin \
	-I$(top_srcdir)/src/utilities/Core \
	-I$(top_srcdir)/src/utilities/Solvers \
	-I$(top_srcdir)/../TSFCoreUtils/src \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
	-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
	-I$(top_builddir)/../TSFCoreUtils/src \
    -I$(top_builddir)/../teuchos/src \
    -I$(top_srcdir)/../teuchos \
    -I$(top_srcdir)/../teuchos/src

TSFCoreNonlinNP4DOptTest_LDADD = \
	$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LIBS) $(FLIBS)


