#@HEADER
# ************************************************************************
# 
#              Meros: Segregated Preconditioning Package
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

EXAMPLE = saddle_kayloghin

noinst_PROGRAMS = saddle_kayloghin

saddle_kayloghin_SOURCES = $(srcdir)/saddle_kayloghin.cpp

saddle_kayloghin_DEPENDENCIES = $(top_builddir)/../TSF/src/libTSF.a

saddle_kayloghin_CXXFLAGS =                       \
-I$(top_builddir)/src                             \
 -I$(top_srcdir)/src                              \
 -I$(top_builddir)/../TSF/src/Support             \
 -I$(top_srcdir)/../TSF/src/Support               \
 -I$(top_srcdir)/../TSF/src/Vectors               \
 -I$(top_srcdir)/../TSF/src/VectorSpaces          \
 -I$(top_srcdir)/../TSF/src/LinearOperators       \
 -I$(top_srcdir)/../TSF/src/Preconditioners       \
 -I$(top_srcdir)/../TSF/src/NonlinearOperators    \
 -I$(top_srcdir)/../TSF/src/LinearSolvers         \
 -I$(top_srcdir)/../TSF/src/NonlinearSolvers      \
 -I$(top_srcdir)/../TSF/src/LinearProblem         \
 -I$(top_srcdir)/../TSF/src/Readers               \
 -I$(top_srcdir)/../TSF/src/Writers               \
 -I$(top_srcdir)/../TSF/src/XML                   \
 -I$(top_srcdir)/../TSF/src/Implementations       \
 -I$(top_builddir)/../ifpack/src                  \
 -I$(top_srcdir)/../ifpack/src                    \
 -I$(top_builddir)/../ml/src                      \
 -I$(top_srcdir)/../ml/src/Coarsen                \
 -I$(top_srcdir)/../ml/src/Comm                   \
 -I$(top_srcdir)/../ml/src/FEGrid                 \
 -I$(top_srcdir)/../ml/src/Include                \
 -I$(top_srcdir)/../ml/src/Krylov                 \
 -I$(top_srcdir)/../ml/src/Main                   \
 -I$(top_srcdir)/../ml/src/MLAPI                  \
 -I$(top_srcdir)/../ml/src/Mrtr                   \
 -I$(top_srcdir)/../ml/src/NonlinML               \
 -I$(top_srcdir)/../ml/src/Operator               \
 -I$(top_srcdir)/../ml/src/Runtime                \
 -I$(top_srcdir)/../ml/src/Smoother               \
 -I$(top_srcdir)/../ml/src/Utils                  \
 -I$(top_builddir)/../epetraext/src               \
 -I$(top_srcdir)/../epetraext/src                 \
 -I$(top_srcdir)/../epetraext/src/block           \
 -I$(top_srcdir)/../epetraext/src/btf             \
 -I$(top_srcdir)/../epetraext/src/coloring        \
 -I$(top_srcdir)/../epetraext/src/distdir         \
 -I$(top_srcdir)/../epetraext/src/inout           \
 -I$(top_srcdir)/../epetraext/src/matlab          \
 -I$(top_srcdir)/../epetraext/src/matlab          \
 -I$(top_srcdir)/../epetraext/src/model_evaluator \
 -I$(top_srcdir)/../epetraext/src/transform       \
 -I$(top_srcdir)/../epetraext/src/zoltan          \
 -I$(top_builddir)/../epetra/src                  \
 -I$(top_srcdir)/../epetra/src                    \
 -I$(top_builddir)/../aztecoo/src                 \
 -I$(top_srcdir)/../aztecoo/src                   \
 -I$(top_builddir)/../teuchos/src                 \
 -I$(top_srcdir)/../teuchos/src                   \
 -I$(top_builddir)/../triutils/src                \
 -I$(top_srcdir)/../triutils/src  

saddle_kayloghin_LDADD =                          \
 $(top_builddir)/../meros/src/libmeros.a          \
 $(top_builddir)/../TSF/src/libTSF.a 	          \
 $(top_builddir)/../ifpack/src/libifpack.a        \
 $(top_builddir)/../ml/src/libml.a                \
 $(top_builddir)/../epetraext/src/libepetraext.a  \
 $(top_builddir)/../epetra/src/libepetra.a        \
 $(top_builddir)/../aztecoo/src/libaztecoo.a      \
 $(top_builddir)/../teuchos/src/libteuchos.a      \
 $(top_builddir)/../triutils/src/libtriutils.a    \
 $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
