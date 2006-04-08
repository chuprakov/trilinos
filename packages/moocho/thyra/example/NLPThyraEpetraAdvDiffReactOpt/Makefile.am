# @HEADER
# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#                  Copyright (2003) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

noinst_PROGRAMS = NLPThyraEpetraAdvDiffReactOpt

NLPThyraEpetraAdvDiffReactOpt_SOURCES = \
  $(srcdir)/GLpApp_AdvDiffReactOptModel.cpp\
  $(srcdir)/GLpApp_GLpYUEpetraDataPool.cpp\
  $(srcdir)/GenSQP_YUEpetraVector.cpp\
  $(srcdir)/assemble.cpp\
  $(srcdir)/assemble_bdry.cpp\
  $(srcdir)/lassembly.cpp\
  $(srcdir)/meshreader.cpp\
  $(srcdir)/nonlinhessvec.cpp\
  $(srcdir)/nonlinjac.cpp\
  $(srcdir)/nonlinvec.cpp\
  $(srcdir)/print2MATLAB.cpp\
  $(srcdir)/quadrature.cpp\
  $(srcdir)/usr_par.cpp\
  $(srcdir)/rect2DMeshGenerator.cpp\
  $(srcdir)/NLPThyraEpetraAdvDiffReactOptMain.cpp

#NLPThyraEpetraAdvDiffReactOpt_DEPENDENCIES = $(common_dependencies) input_mesh_files
NLPThyraEpetraAdvDiffReactOpt_DEPENDENCIES = $(common_dependencies) square.1.000
NLPThyraEpetraAdvDiffReactOpt_LDADD = $(common_ldadd)


square.1.000: $(srcdir)/data/square/square.1.000
	cp $(srcdir)/data/square/square.1.000 square.1.000
	cp $(srcdir)/data/square/square.2.000 square.2.000
	cp $(srcdir)/data/square/square.9x9.2.000 square.9x9.2.000
	cp $(srcdir)/data/square/square.9x9.2.001 square.9x9.2.001

include $(top_srcdir)/MakefileExeStub.mak

CPPFLAGS += -I$(top_srcdir)/../thyra/example/operator_vector
