#@HEADER
#@HEADER

EXEEXT = .exe

EXAMPLE = saddle_kayloghin

noinst_PROGRAMS = saddle_kayloghin

saddle_kayloghin_SOURCES = $(srcdir)/saddle_kayloghin.cpp

saddle_kayloghin_DEPENDENCIES = $(top_builddir)/../TSF/src/libTSF.a

saddle_kayloghin_CXXFLAGS =           -I$(top_srcdir)/src  \
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

saddle_kayloghin_LDADD =  $(top_builddir)/src/libmeros.a \
 $(top_builddir)/../TSF/src/libTSF.a 	         \
 $(top_builddir)/../ifpack/src/libifpack.a       \
 $(top_builddir)/../epetra/src/libepetra.a       \
 $(top_builddir)/../triutils/src/libtriutils.a   \
 $(top_builddir)/../ml/Obj/libml.a               \
 $(top_builddir)/../aztecoo/src/libaztecoo.a  $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
