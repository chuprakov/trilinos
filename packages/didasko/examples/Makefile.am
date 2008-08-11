#@HEADER
# ************************************************************************
# 
#               Disasko : a Tutorial for Trilinos
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

# The following line helps the test harness recover from build errors.

all-local:

if HAVE_DIDASKO_AMESOS
AMESOS_SUBDIR = amesos
endif

if HAVE_DIDASKO_ANASAZI
ANASAZI_SUBDIR = anasazi
endif

if HAVE_DIDASKO_AZTECOO
AZTECOO_SUBDIR = aztecoo
endif

if HAVE_DIDASKO_BELOS
BELOS_SUBDIR = belos
endif

if HAVE_DIDASKO_EPETRA
EPETRA_SUBDIR = epetra
endif

if HAVE_DIDASKO_EPETRAEXT
EPETRAEXT_SUBDIR = epetraext
endif

if HAVE_DIDASKO_IFPACK
IFPACK_SUBDIR = ifpack
endif

if HAVE_DIDASKO_ML
ML_SUBDIR = ml
endif

if HAVE_DIDASKO_NOX
NOX_SUBDIR = nox
endif

if HAVE_DIDASKO_TEUCHOS
TEUCHOS_SUBDIR = teuchos
endif

if HAVE_DIDASKO_TRIUTILS
TRIUTILS_SUBDIR = triutils
endif

if HAVE_DIDASKO_TPETRA
TPETRA_SUBDIR = tpetra
endif


SUBDIRS = $(AMESOS_SUBDIR) $(ANASAZI_SUBDIR) $(AZTECOO_SUBDIR) \
	  $(BELOS_SUBDIR) $(EPETRA_SUBDIR) $(EPETRAEXT_SUBDIR) \
	  $(IFPACK_SUBDIR)  $(ML_SUBDIR) $(NOX_SUBDIR) $(TEUCHOS_SUBDIR) \
	  $(TRIUTILS_SUBDIR) $(TPETRA_SUBDIR) HBMatrices
