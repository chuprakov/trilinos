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

%module(package="PyTrilinos.LOCA.Pitchfork") MinimallyAugmented

%{
// LOCA includes
#include "LOCA.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP handling
%teuchos_rcp(LOCA::Pitchfork::MinimallyAugmented::AbstractGroup)

// Base class support
%pythoncode
%{
import os
import sys
parentDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.normpath(os.path.join(parentDir,".."))
if not parentDir in sys.path: sys.path.append(parentDir)
%}
// %import "NOX.Abstract.i"
%import "LOCA.TurningPoint.MinimallyAugmented_RelPath.i"
%import "LOCA.Pitchfork.MooreSpence.i"

// LOCA::Pitchfork::MinimallyAugmented AbstractGroup class
%include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
