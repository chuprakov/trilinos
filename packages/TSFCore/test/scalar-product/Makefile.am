#@HEADER
# ************************************************************************
# 
#               TSFCore: Trilinos Solver Framework Core
#                 Copyright (2004) Sandia Corporation
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

EXEEXT = .exe

if BUILD_TEUCHOS
SCALAR_PRODUCT_TEST = scalar_product
endif

noinst_PROGRAMS = $(SCALAR_PRODUCT_TEST)

scalar_product_SOURCES = $(srcdir)/cxx_main.cpp

scalar_product_DEPENDENCIES = \
  $(top_builddir)/src/libTSFCore.a

scalar_product_LDFLAGS =

scalar_product_CXXFLAGS = \
  -I$(top_builddir)/src \
  -I$(top_srcdir)/src \
  -I$(top_srcdir)/src/interfaces/Core \
  -I$(top_srcdir)/src/utilities/Core \
  -I$(top_builddir)/../TSFCoreUtils/src \
  -I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
  -I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
  -I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
  -I$(top_srcdir)/../TSFCoreUtils/src \
  -I$(top_builddir)/../teuchos/src \
  -I$(top_srcdir)/../teuchos/src

scalar_product_LDADD = \
  $(top_builddir)/src/libTSFCore.a \
  $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
  $(top_builddir)/../teuchos/src/libteuchos.a \
  $(BLAS_LIBS) $(FLIBS)
