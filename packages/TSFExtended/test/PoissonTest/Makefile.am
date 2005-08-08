#@HEADER
# ************************************************************************
# 
#           TSFExtended: Trilinos Solver Framework Extended
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

include $(top_builddir)/Makefile.export.TSFExtended


EXEEXT = .exe

noinst_PROGRAMS = Poisson_test AztecPoisson ParameterListAztecPoisson



AztecPoisson_SOURCES = $(srcdir)/AztecPoisson.cpp 


AztecPoisson_DEPENDENCIES = $(TSFEXTENDED_LIBS)

AztecPoisson_CXXFLAGS = $(CXXFLAGS) $(EXPORT_INCLUDES)


AztecPoisson_LDADD = $(LDFLAGS) $(TSFEXTENDED_LIBS)





ParameterListAztecPoisson_SOURCES = $(srcdir)/ParameterListAztecPoisson.cpp 


ParameterListAztecPoisson_DEPENDENCIES = $(TSFEXTENDED_LIBS)

ParameterListAztecPoisson_CXXFLAGS = $(CXXFLAGS) $(EXPORT_INCLUDES)

ParameterListAztecPoisson_LDADD = $(LDFLAGS) $(TSFEXTENDED_LIBS)





Poisson_test_SOURCES = $(srcdir)/cxx_main.cpp 

Poisson_test_DEPENDENCIES = $(TSFEXTENDED_LIBS)

Poisson_test_CXXFLAGS = $(CXXFLAGS) $(EXPORT_INCLUDES)

Poisson_test_LDADD = $(LDFLAGS) $(TSFEXTENDED_LIBS)
