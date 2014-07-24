// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
// Teuchos include
#include "PyTrilinos_Teuchos_Util.h"

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "Epetra_NumPyVector.h"
%}
#endif

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::Hopf::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::ExtendedGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::SolverStrategy)
%teuchos_rcp(LOCA::Hopf::MooreSpence::SolverFactory)
%teuchos_rcp(LOCA::Hopf::MooreSpence::SalingerBordering)

// Base class support
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}
%import "LOCA.Extended_RelPath.i"
%import "LOCA.MultiContinuation_RelPath.i"
%import "LOCA.TimeDependent_RelPath.i"
%import "LOCA.TurningPoint.MooreSpence_RelPath.i"

// LOCA::Hopf::MooreSpence AbstractGroup class
%include "LOCA_Hopf_MooreSpence_AbstractGroup.H"

// LOCA::Hopf::MooreSpence ExtendedGroup class
%include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"

// LOCA::Hopf::MooreSpence ExtendedMultiVector class
%include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"

// LOCA::Hopf::MooreSpence ExtendedVector class
%ignore LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency;
%ignore LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam;
%include "LOCA_Hopf_MooreSpence_ExtendedVector.H"

// LOCA::Hopf::MooreSpence FiniteDifferenceGroup class
%include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"

// LOCA::Hopf::MooreSpence SolverFactory class
%include "LOCA_Hopf_MooreSpence_SolverFactory.H"

// LOCA::Hopf::MooreSpence SolverStrategy class
%include "LOCA_Hopf_MooreSpence_SolverStrategy.H"

// LOCA::Hopf::MooreSpence SalingerBordering class
%include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
