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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

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
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup)

// Base class support
%import "LOCA.TurningPoint.MooreSpence.i"

// LOCA::TurningPoint::MinimallyAugmented AbtractGroup class
%feature("director") LOCA::TurningPoint::MinimallyAugmented;
%include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"

// LOCA::TurningPoint::MinimallyAugmented FinitDifferenceGroup class
%feature("director") LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup;
%include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
