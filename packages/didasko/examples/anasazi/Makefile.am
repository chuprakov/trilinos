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

include $(top_builddir)/Makefile.export

EXAMPLE = ex1 ex2 ex3 ex4
EXEEXT = .exe

noinst_PROGRAMS = $(EXAMPLE)

ex1_SOURCES = $(srcdir)/ex1.cpp
ex1_CXXFLAGS = $(DIDASKO_INCLUDES) 
ex1_LDADD = $(DIDASKO_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex2_SOURCES = $(srcdir)/ex2.cpp
ex2_CXXFLAGS = $(DIDASKO_INCLUDES) 
ex2_LDADD = $(DIDASKO_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex3_SOURCES = $(srcdir)/ex3.cpp
ex3_CXXFLAGS = $(DIDASKO_INCLUDES) 
ex3_LDADD = $(DIDASKO_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)

ex4_SOURCES = $(srcdir)/ex4.cpp
ex4_CXXFLAGS = $(DIDASKO_INCLUDES) 
ex4_LDADD = $(DIDASKO_LIBS) \
            $(LIBS)         \
            $(LAPACK_LIBS)  \
            $(BLAS_LIBS)    \
            $(FLIBS)
