
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

#JW added all paths beginning with -I$(top_builddir) so that TSF could be
#configured and built in a build tree separate from the source tree.
AM_CPPFLAGS = $(TSF_INCLUDES)	-I$(top_srcdir)/../aztecoo/src -I$(top_srcdir)/../ifpack/src -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src -I$(top_builddir)/../ifpack/src -I$(top_builddir)/../aztecoo/src

#JW changed top_srcdir to top_builddir in each line below.  If source and 
#build directories are not the same, the libraries will be located in the
#build directory.
TSF_LIBS = -L$(top_builddir)/src \
	-L$(top_builddir)/../aztecoo/src \
	-L$(top_builddir)/../ifpack/src \
	-L$(top_builddir)/../epetra/src \
	-L$(top_builddir)/../triutils/src \
	-ltsf -laztecoo -lifpack -lepetra -ltriutils -ly12m

