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

#
# Actual programs
#

SILLY_CG_SOLVE_SERIAL = sillyCgSolve_serial
SILLY_POWER_METHOD_SERIAL = sillyPowerMethod_serial
if HAVE_TSFCORE_EPETRA
SILLY_CG_SOLVE_EPETRA = sillyCgSolve_epetra
SILLY_POWER_METHOD_EPETRA = sillyPowerMethod_epetra
else
SILLY_CG_SOLVE_EPETRA = 
SILLY_POWER_METHOD_EPETRA = 
endif

noinst_PROGRAMS = $(SILLY_CG_SOLVE_SERIAL) $(SILLY_CG_SOLVE_EPETRA) $(SILLY_POWER_METHOD_SERIAL) $(SILLY_POWER_METHOD_EPETRA)

sillyCgSolve_serial_SOURCES = $(srcdir)/sillyCgSolve_serial.cpp
sillyCgSolve_serial_DEPENDENCIES = $(common_dependencies)
sillyCgSolve_serial_CXXFLAGS = $(common_cxxflags)
sillyCgSolve_serial_LDADD = $(common_ldadd)

sillyPowerMethod_serial_SOURCES = $(srcdir)/sillyPowerMethod_serial.cpp
sillyPowerMethod_serial_DEPENDENCIES = $(common_dependencies)
sillyPowerMethod_serial_CXXFLAGS = $(common_cxxflags)
sillyPowerMethod_serial_LDADD = $(common_ldadd)

sillyPowerMethod_epetra_SOURCES = $(srcdir)/sillyPowerMethod_epetra.cpp \
                                  $(srcdir)/createTridiagEpetraLinearOp.cpp
sillyPowerMethod_epetra_DEPENDENCIES = $(common_dependencies)
sillyPowerMethod_epetra_CXXFLAGS = $(common_cxxflags)
sillyPowerMethod_epetra_LDADD = $(common_ldadd)

sillyCgSolve_epetra_SOURCES = $(srcdir)/sillyCgSolve_epetra.cpp \
                              $(srcdir)/createTridiagEpetraLinearOp.cpp
sillyCgSolve_epetra_DEPENDENCIES = $(common_dependencies)
sillyCgSolve_epetra_CXXFLAGS = $(common_cxxflags)
sillyCgSolve_epetra_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
	$(top_builddir)/src/libTSFCore.a \
	$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a

common_cxxflags = \
	-I$(top_builddir)/../TSFCore/src \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src/interfaces/Core \
	-I$(top_srcdir)/src/utilities/Core \
	-I$(top_srcdir)/src/adapters/serial \
	-I$(top_srcdir)/src/adapters/mpi-base \
	-I$(top_srcdir)/src/adapters/epetra \
	-I$(top_builddir)/../TSFCoreUtils/src \
	-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
	-I$(top_srcdir)/../TSFCoreUtils/src \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
	-I$(top_builddir)/../teuchos/src \
	-I$(top_srcdir)/../teuchos \
	-I$(top_srcdir)/../teuchos/src

if HAVE_TSFCORE_EPETRA
common_cxxflags += \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../epetra/src
endif

common_ldadd = \
	$(top_builddir)/src/libTSFCore.a \
	$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

if HAVE_TSFCORE_EPETRA
common_ldadd += \
	$(top_builddir)/../epetra/src/libepetra.a
endif

common_ldadd += \
	$(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)










