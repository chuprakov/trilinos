
EXEEXT = .exe

noinst_PROGRAMS = cgTest
cgTest_SOURCES = Tests/cgTest.cpp

cgTest_LDADD = $(TSF_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

cgTest_DEPENDENCIES = $(top_builddir)/src/libtsf.a

cgTest_CXXFLAGS = $(SUNDANCE_INCLUDES)

%.e : %.o
	$(CXXLINK) -o $@ $(LDFLAGS) $< -L$(top_builddir)/src $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

%.o : %.cpp
	$(CXXCOMPILE) -c -o $@ $< 





TSF_INCLUDES = -I$(top_srcdir)/src/Support -I$(top_srcdir)/src/Vectors \
	-I$(top_srcdir)/src/VectorSpaces -I$(top_srcdir)/src/LinearOperators \
	-I$(top_srcdir)/src/Preconditioners -I$(top_srcdir)/src/NonlinearOperators \
	-I$(top_srcdir)/src/LinearSolvers -I$(top_srcdir)/src/NonlinearSolvers \
	-I$(top_srcdir)/src/LinearProblem -I$(top_srcdir)/src/Readers \
	-I$(top_srcdir)/src/Writers -I$(top_srcdir)/src/XML \
	-I$(top_srcdir)/src/Implementations

AM_CPPFLAGS = $(TSF_INCLUDES)	-I$(top_srcdir)/../aztecoo/src -I$(top_srcdir)/../ifpack/src -I$(top_srcdir)/../epetra/src 


TSF_LIBS = -L$(top_srcdir)/src \
	-L$(top_srcdir)/../aztecoo/src \
	-L$(top_srcdir)/../ifpack/src \
	-L$(top_srcdir)/../epetra/src \
	-L$(top_srcdir)/../triutils/src \
	-ltsf -laztecoo -lifpack -lepetra -ltriutils -ly12m

