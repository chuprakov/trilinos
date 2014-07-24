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

%define %loca_abstract_docstring
"
PyTrilinos.LOCA.Abstract is the python interface to namespace Abstract
of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Abstract is to provide abstract continuation
problem base classes.  The python version of LOCA.Abstract supports
the following classes:

    * Group                - 
    * TransposeSolveGroup  - 
    * Iterator             - 
    * Factory              - 

Any other notes about the package as a whole. . . .
"
%enddef

%module(package      = "PyTrilinos.LOCA",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_abstract_docstring) Abstract

%{
// PyTrilinos includes
#include "PyTrilinos_PythonException.h"
#include "PyTrilinos_Teuchos_Util.h"

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"
// #include "LOCA_MultiContinuation_AbstractStrategy.H"
// #include "LOCA_MultiContinuation_ArcLengthConstraint.H"
// #include "LOCA_MultiContinuation_ArcLengthGroup.H"
// #include "LOCA_MultiContinuation_ExtendedGroup.H"
// #include "LOCA_MultiContinuation_CompositeConstraint.H"
// #include "LOCA_MultiContinuation_CompositeConstraintMVDX.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Director exception handling
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(PyTrilinos::PythonException & e)
  {
    e.restore();
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore *::operator=;
%ignore *::operator[];

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP handling
%teuchos_rcp(LOCA::Abstract::Group)
%teuchos_rcp(LOCA::Abstract::TransposeSolveGroup)
%teuchos_rcp(LOCA::Abstract::Iterator)
%teuchos_rcp(LOCA::Abstract::Factory)

// Import SWIG interface files to provide information about base
// classes
%import "LOCA.Homotopy.i"
%import "LOCA.PhaseTransition.i"
%import "LOCA.TurningPoint.MinimallyAugmented.i"
%import "LOCA.Hopf.MinimallyAugmented.i"
%import "LOCA.Pitchfork.MinimallyAugmented.i"

// LOCA::Abstract Group class
%include "LOCA_Abstract_Group.H"

// LOCA::Abstract TransposeSolveGroup class
%include "LOCA_Abstract_TransposeSolveGroup.H"

// LOCA::Abstract Iterator class
%include "LOCA_Abstract_Iterator.H"

// LOCA::Abstract Factory class
%include "LOCA_Abstract_Factory.H"
