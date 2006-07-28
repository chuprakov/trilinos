include $(top_builddir)/Makefile.export

EXAMPLE = ex_Comm \
	  ex_Map \
	  ex_Vector \
	  ex_MultiVector \
	  ex_CrsMatrix \
	  ex_FECrsMatrix \
	  ex_VbrMatrix
EXEEXT = .exe

noinst_PROGRAMS = $(EXAMPLE)

ex_Comm_SOURCES = $(srcdir)/ex_Comm.cpp
ex_Comm_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_Comm_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex_Map_SOURCES = $(srcdir)/ex_Map.cpp
ex_Map_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_Map_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex_Vector_SOURCES = $(srcdir)/ex_Vector.cpp
ex_Vector_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_Vector_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex_MultiVector_SOURCES = $(srcdir)/ex_MultiVector.cpp
ex_MultiVector_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_MultiVector_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)
            $(FLIBS)

ex_CrsMatrix_SOURCES = $(srcdir)/ex_CrsMatrix.cpp
ex_CrsMatrix_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_CrsMatrix_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex_FECrsMatrix_SOURCES = $(srcdir)/ex_FECrsMatrix.cpp
ex_FECrsMatrix_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_FECrsMatrix_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex_VbrMatrix_SOURCES = $(srcdir)/ex_VbrMatrix.cpp
ex_VbrMatrix_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_VbrMatrix_LDADD = $(WEBTRILINOS_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)
