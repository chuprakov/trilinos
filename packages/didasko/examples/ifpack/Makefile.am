#@HEADER
# ************************************************************************
# 
#               Didasko : a Tutorial for Trilinos
#                 Copyright (2001) Sandia Corporation
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

EXAMPLE = ex1 ex2
EXEEXT = .exe

DIDASKO_INCLUDES = -I$(top_srcdir)/src

if HAVE_DIDASKO_EPETRA
include $(top_builddir)/../epetra/Makefile.export
else
EPETRA_INCLUDE = 
EPETRA_LIBS =
endif

if HAVE_DIDASKO_AZTECOO
AZTECOO_INCLUDES = -I$(top_builddir)/../aztecoo/src \
                  -I$(top_srcdir)/../aztecoo/src
AZTECOO_LIBS = $(top_builddir)/../aztecoo/src/libaztecoo.a
else
AZTECOO_INCLUDES = 
AZTECOO_LIBS =
endif

if HAVE_DIDASKO_IFPACK
IFPACK_INCLUDES = -I$(top_builddir)/../ifpack/src \
                  -I$(top_srcdir)/../ifpack/src
IFPACK_LIBS = $(top_builddir)/../ifpack/src/libifpack.a
else
IFPACK_INCLUDES = 
IFPACK_LIBS =
endif

if HAVE_DIDASKO_TEUCHOS
TEUCHOS_INCLUDES = -I$(top_builddir)/../teuchos/src \
                  -I$(top_srcdir)/../teuchos/src
TEUCHOS_LIBS = $(top_builddir)/../teuchos/src/libteuchos.a
else
TEUCHOS_INCLUDES = 
TEUCHOS_LIBS =
endif

if HAVE_DIDASKO_TRIUTILS
TRIUTILS_INCLUDES = -I$(top_builddir)/../triutils/src \
                  -I$(top_srcdir)/../triutils/src
TRIUTILS_LIBS = $(top_builddir)/../triutils/src/libtriutils.a
else
TRIUTILS_INCLUDES = 
TRIUTILS_LIBS =
endif

noinst_PROGRAMS = $(EXAMPLE)

ex1_SOURCES = $(srcdir)/ex1.cpp
ex1_DEPENDENCIES = $(EPETRA_LIBS) \
		   $(AZTECOO_LIBS) \
		   $(IFPACK_LIBS) \
		   $(TRIUTILS_LIBS) \
		   $(TEUCHOS_LIBS)
ex1_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(EPETRA_INCLUDE) \
	       $(AZTECOO_INCLUDES) \
	       $(TRIUTILS_INCLUDES) \
	       $(IFPACK_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex1_LDADD = $(IFPACK_LIBS) \
	    $(TRIUTILS_LIBS) \
	    $(AZTECOO_LIBS) \
	    $(EPETRA_LIBS) \
	    $(TEUCHOS_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)

ex2_SOURCES = $(srcdir)/ex2.cpp
ex2_DEPENDENCIES = $(EPETRA_LIBS) \
		   $(AZTECOO_LIBS) \
		   $(IFPACK_LIBS) \
		   $(TRIUTILS_LIBS) \
		   $(TEUCHOS_LIBS)
ex2_CXXFLAGS = $(DIDASKO_INCLUDES) \
	       $(EPETRA_INCLUDE) \
	       $(AZTECOO_INCLUDES) \
	       $(TRIUTILS_INCLUDES) \
	       $(IFPACK_INCLUDES) \
	       $(TEUCHOS_INCLUDES)
ex2_LDADD = $(IFPACK_LIBS) \
	    $(TRIUTILS_LIBS) \
	    $(AZTECOO_LIBS) \
	    $(EPETRA_LIBS) \
	    $(TEUCHOS_LIBS) \
	    $(LAPACK_LIBS) $(LIBS) $(BLAS_LIBS) $(FLIBS)
