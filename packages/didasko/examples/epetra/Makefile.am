#@HEADER
# ************************************************************************
# 
#               Didasko : a Tutorial for Trilinos
#                 Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXAMPLE = ex1  ex2  ex3  ex4  ex5  ex6  ex7  ex8  ex9  ex10 \
	  ex11 ex12 ex13 ex14 ex15 ex16 ex17 ex18 ex19 ex20 \
	  ex21 ex22 ex23 ex24 ex25
EXEEXT = .exe

DIDASKO_INCLUDES = -I$(top_srcdir)/src

if HAVE_DIDASKO_EPETRA
include $(top_builddir)/../epetra/Makefile.export
else
EPETRA_INCLUDE = 
EPETRA_LIBS = =
endif

# Triutils is used by example 25 only
if HAVE_DIDASKO_TRIUTILS
TRIUTILS_INCLUDES = -I$(top_srcdir)/../triutils/src \
                  -I$(top_builddir)/../triutils/src
TRIUTILS_LIBS = $(top_builddir)/../triutils/src/libtriutils.a
else
TRIUTILS_INCLUDES = 
TRIUTILS_LIBS = =
endif
noinst_PROGRAMS = $(EXAMPLE)

ex1_SOURCES = $(srcdir)/ex1.cpp
ex1_DEPENDENCIES = $(EPETRA_LIBS)
ex1_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex1_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex2_SOURCES = $(srcdir)/ex2.cpp
ex2_DEPENDENCIES = $(EPETRA_LIBS)
ex2_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex2_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex3_SOURCES = $(srcdir)/ex3.cpp
ex3_DEPENDENCIES = $(EPETRA_LIBS)
ex3_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex3_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex4_SOURCES = $(srcdir)/ex4.cpp
ex4_DEPENDENCIES = $(EPETRA_LIBS)
ex4_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex4_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex5_SOURCES = $(srcdir)/ex5.cpp
ex5_DEPENDENCIES = $(EPETRA_LIBS)
ex5_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex5_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex6_SOURCES = $(srcdir)/ex6.cpp
ex6_DEPENDENCIES = $(EPETRA_LIBS)
ex6_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex6_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex7_SOURCES = $(srcdir)/ex7.cpp
ex7_DEPENDENCIES = $(EPETRA_LIBS)
ex7_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex7_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex8_SOURCES = $(srcdir)/ex8.cpp
ex8_DEPENDENCIES = $(EPETRA_LIBS)
ex8_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex8_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex9_SOURCES = $(srcdir)/ex9.cpp
ex9_DEPENDENCIES = $(EPETRA_LIBS)
ex9_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex9_LDADD = $(EPETRA_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex10_SOURCES = $(srcdir)/ex10.cpp
ex10_DEPENDENCIES = $(EPETRA_LIBS)
ex10_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex10_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex11_SOURCES = $(srcdir)/ex11.cpp
ex11_DEPENDENCIES = $(EPETRA_LIBS)
ex11_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex11_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex12_SOURCES = $(srcdir)/ex12.cpp
ex12_DEPENDENCIES = $(EPETRA_LIBS)
ex12_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex12_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex13_SOURCES = $(srcdir)/ex13.cpp
ex13_DEPENDENCIES = $(EPETRA_LIBS)
ex13_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex13_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex14_SOURCES = $(srcdir)/ex14.cpp
ex14_DEPENDENCIES = $(EPETRA_LIBS)
ex14_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex14_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex15_SOURCES = $(srcdir)/ex15.cpp
ex15_DEPENDENCIES = $(EPETRA_LIBS)
ex15_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex15_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex16_SOURCES = $(srcdir)/ex16.cpp
ex16_DEPENDENCIES = $(EPETRA_LIBS)
ex16_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex16_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex17_SOURCES = $(srcdir)/ex17.cpp
ex17_DEPENDENCIES = $(EPETRA_LIBS)
ex17_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex17_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex18_SOURCES = $(srcdir)/ex18.cpp
ex18_DEPENDENCIES = $(EPETRA_LIBS)
ex18_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex18_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex19_SOURCES = $(srcdir)/ex19.cpp
ex19_DEPENDENCIES = $(EPETRA_LIBS)
ex19_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex19_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex20_SOURCES = $(srcdir)/ex20.cpp
ex20_DEPENDENCIES = $(EPETRA_LIBS)
ex20_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex20_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex21_SOURCES = $(srcdir)/ex21.cpp
ex21_DEPENDENCIES = $(EPETRA_LIBS)
ex21_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex21_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex22_SOURCES = $(srcdir)/ex22.cpp
ex22_DEPENDENCIES = $(EPETRA_LIBS)
ex22_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex22_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex23_SOURCES = $(srcdir)/ex23.cpp
ex23_DEPENDENCIES = $(EPETRA_LIBS)
ex23_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex23_LDADD = $(EPETRA_LIBS) \
	      $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex24_SOURCES = $(srcdir)/ex24.cpp
ex24_DEPENDENCIES = $(EPETRA_LIBS)
ex24_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE)
ex24_LDADD = $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex25_SOURCES = $(srcdir)/ex25.cpp
ex25_DEPENDENCIES = $(EPETRA_LIBS) $(TRIUTILS_LIBS)
ex25_CXXFLAGS = $(DIDASKO_INCLUDES) $(EPETRA_INCLUDE) $(TRIUTILS_INCLUDES)
ex25_LDADD = $(TRIUTILS_LIBS) \
	     $(EPETRA_LIBS) \
	     $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)
