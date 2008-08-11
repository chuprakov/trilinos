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

%define %nox_abstract_docstring
"
PyTrilinos.NOX.Abstract is the python interface to namespace Abstract
of the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.Abstract is to provide base classes from which
concrete NOX interfaces can be derived.  Currently, the only concrete
implementation is for Epetra, in the NOX.Epetra module.

NOX.Abstract provides the following user-level classes:

    * Group            - Class defining a collection of objects needed by NOX
    * PrePostOperator  - Pre- and post-iteration operators
    * MultiVector      - Multivector class
    * Vector           - Vector class
"
%enddef

%module(package      = "PyTrilinos.NOX",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_abstract_docstring) Abstract

%{
// Teuchos includes
#include "Teuchos_PythonParameter.h"

// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Abstract_Vector.H"

// Local includes
#include "NumPyImporter.h"
%}

// Standard exception handling
%include "exception.i"

// Include NOX documentation
%include "NOX_dox.i"

// General ignore directives
%ignore *::operator=;
%ignore *::operator[];

// Trilinos module imports
%import "Teuchos.i"

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unkown C++ exception");
  }
}

// Support for Teuchos::RCPs
%teuchos_rcp_typemaps(NOX::Abstract::Group)

////////////////////////////////
// NOX_Abstract_Group support //
////////////////////////////////
%include "NOX_Abstract_Group.H"

//////////////////////////////////////////
// NOX_Abstract_PrePostOperator support //
//////////////////////////////////////////
%feature("director") NOX::Abstract::PrePostOperator;
%include "NOX_Abstract_PrePostOperator.H"

//////////////////////////////////////
// NOX_Abstract_MultiVector support //
//////////////////////////////////////
%ignore NOX::Abstract::MultiVector::clone(int) const;
%rename(_print) NOX::Abstract::MultiVector::print;
%include "NOX_Abstract_MultiVector.H"

/////////////////////////////////
// NOX_Abstract_Vector support //
/////////////////////////////////
%rename(_print) NOX::Abstract::Vector::print;
%include "NOX_Abstract_Vector.H"

// Turn off the exception handling
%exception;
