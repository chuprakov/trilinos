#@HEADER
#@HEADER

EXEEXT = .exe

EXAMPLE = saddle_solve

noinst_PROGRAMS = saddle_solve

saddle_solve_SOURCES = $(srcdir)/saddle_solve.cpp

saddle_solve_DEPENDENCIES = $(top_builddir)/../TSF/src/libTSF.a

saddle_solve_CXXFLAGS =           -I$(top_srcdir)/src  \
          -I $(top_srcdir)/../TSF/src/LinearOperators  \
          -I$(top_srcdir)/../TSF/src/Implementations  \
          -I$(top_srcdir)/../TSF/src/LinearProblem    \
          -I$(top_srcdir)/../TSF/src/Vectors          \
          -I$(top_srcdir)/../TSF/src/Support          \
          -I$(top_srcdir)/../TSF/src/Preconditioners  \
          -I$(top_srcdir)/../TSF/src/VectorSpaces     \
          -I$(top_srcdir)/../TSF/src/LinearSolvers    \
          -I$(top_srcdir)/../aztecoo/src              \
          -I$(top_srcdir)/../epetra/src               \
          -I$(top_srcdir)/../ml/Obj            \
          -I$(top_builddir)/src  \
          -I $(top_builddir)/../TSF/src/LinearOperators  \
          -I$(top_builddir)/../TSF/src/Implementations  \
          -I$(top_builddir)/../TSF/src/LinearProblem    \
          -I$(top_builddir)/../TSF/src/Vectors          \
          -I$(top_builddir)/../TSF/src/Support          \
          -I$(top_builddir)/../TSF/src/Preconditioners  \
          -I$(top_builddir)/../TSF/src/VectorSpaces     \
          -I$(top_builddir)/../TSF/src/LinearSolvers    \
          -I$(top_builddir)/../aztecoo/src              \
          -I$(top_builddir)/../epetra/src               \
          -I$(top_builddir)/../ml/Obj            

saddle_solve_LDADD =  $(top_builddir)/src/libmeros.a \
 $(top_builddir)/../TSF/src/libTSF.a 	         \
 $(top_builddir)/../ifpack/src/libifpack.a       \
 $(top_builddir)/../epetra/src/libepetra.a       \
 $(top_builddir)/../triutils/src/libtriutils.a   \
 $(top_builddir)/../ml/Obj/libml.a               \
 $(top_builddir)/../aztecoo/src/libaztecoo.a  $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
