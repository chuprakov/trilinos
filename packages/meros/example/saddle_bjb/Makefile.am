#@HEADER
#@HEADER

EXEEXT = .exe

EXAMPLE = saddle_bjb

noinst_PROGRAMS = saddle_bjb

saddle_bjb_SOURCES = $(srcdir)/saddle_bjb.cpp

saddle_bjb_DEPENDENCIES = $(top_builddir)/../TSF/src/libTSF.a

saddle_bjb_CXXFLAGS =           -I$(top_srcdir)/src  \
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
          -I$(top_srcdir)/../ml/src/                  \
          -I$(top_srcdir)/../ml/src/Coarsen           \
          -I$(top_srcdir)/../ml/src/Comm              \
          -I$(top_srcdir)/../ml/src/FEGrid            \
          -I$(top_srcdir)/../ml/src/Include           \
          -I$(top_srcdir)/../ml/src/Krylov            \
          -I$(top_srcdir)/../ml/src/Main              \
          -I$(top_srcdir)/../ml/src/Operator          \
          -I$(top_srcdir)/../ml/src/Smoother          \
          -I$(top_srcdir)/../ml/src/Utils             \
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
          -I$(top_builddir)/../ml/src/                  \
          -I$(top_builddir)/../ml/src/Coarsen           \
          -I$(top_builddir)/../ml/src/Comm              \
          -I$(top_builddir)/../ml/src/FEGrid            \
          -I$(top_builddir)/../ml/src/Include           \
          -I$(top_builddir)/../ml/src/Krylov            \
          -I$(top_builddir)/../ml/src/Main              \
          -I$(top_builddir)/../ml/src/Operator          \
          -I$(top_builddir)/../ml/src/Smoother          \
          -I$(top_builddir)/../ml/src/Utils


saddle_bjb_LDADD =  $(top_builddir)/src/libmeros.a \
 $(top_builddir)/../TSF/src/libTSF.a 	         \
 $(top_builddir)/../ifpack/src/libifpack.a       \
 $(top_builddir)/../epetra/src/libepetra.a       \
 $(top_builddir)/../triutils/src/libtriutils.a   \
 $(top_builddir)/../ml/src/libml.a               \
 $(top_builddir)/../aztecoo/src/libaztecoo.a  $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
