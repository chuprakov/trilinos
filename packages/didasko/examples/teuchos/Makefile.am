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

EXAMPLE = ex1  ex2  ex3  ex4  ex5  ex6  ex7
EXEEXT = .exe

DIDASKO_INCLUDES = -I$(top_srcdir)/src

if HAVE_DIDASKO_TEUCHOS
TEUCHOS_INCLUDES = -I$(top_srcdir)/../teuchos/src \
                 -I$(top_builddir)/../teuchos/src
TEUCHOS_LIBS = $(top_builddir)/../teuchos/src/libteuchos.a
else
TEUCHOS_INCLUDES =
TEUCHOS_LIBS =
endif

noinst_PROGRAMS = $(EXAMPLE)

ex1_SOURCES = $(srcdir)/ex1.cpp
ex1_DEPENDENCIES = $(TEUCHOS_LIBS)
ex1_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex1_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex2_SOURCES = $(srcdir)/ex2.cpp
ex2_DEPENDENCIES = $(TEUCHOS_LIBS)
ex2_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex2_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex3_SOURCES = $(srcdir)/ex3.cpp
ex3_DEPENDENCIES = $(TEUCHOS_LIBS)
ex3_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex3_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex4_SOURCES = $(srcdir)/ex4.cpp
ex4_DEPENDENCIES = $(TEUCHOS_LIBS)
ex4_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex4_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex5_SOURCES = $(srcdir)/ex5.cpp
ex5_DEPENDENCIES = $(TEUCHOS_LIBS)
ex5_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex5_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex6_SOURCES = $(srcdir)/ex6.cpp
ex6_DEPENDENCIES = $(TEUCHOS_LIBS)
ex6_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex6_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex7_SOURCES = $(srcdir)/ex7.cpp
ex7_DEPENDENCIES = $(TEUCHOS_LIBS)
ex7_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex7_LDADD = $(TEUCHOS_LIBS) \
            $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)
