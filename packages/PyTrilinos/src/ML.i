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

%define %ml_docstring
"
PyTrilinos.ML is the python interface to the Trilinos multi-level
preconditioner package ML/MLAPI:

    http://trilinos.sandia.gov/packages/ml

The purpose of ML is to provide multilevel preconditioners to
Trilinos.

ML provides the following user-level classes:

    * BaseObject                - Base class for MLAPI objects
    * CompObject                - FLOP counting base class
    * TimeObject                - Timing base class
    * MultiLevelPreconditioner  - Black-box multilevel preconditioner
    * EpetraBaseOperator        - Base class for interface to Epetra
    * BaseOperator              - Base class for all MLAPI operators
    * Space                     - Defines number of elements and their distribution
    * MultiVector               - MLAPI multivector class
    * Operator                  - MLAPI operator class
    * InverseOperator           - MLAPI inverse operator class
    * PyMatrix                  - Python interface to MLAPI operators

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exMLAPI.py
    * exMLAPI_Simple.py
    * exMLAPI_AztecOO.py
    * exMLAPI_Iterate.py
    * exMLAPI_PyMatrix.py
    * exMLAPI_Smoother.py
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	autodoc   = "1",
	docstring = %ml_docstring) ML

// This is to avoid BaseLinearCombination and derived classes.
// MLAPI_LC is defined in setup.py; the MLAPI code contains some 
// `#ifndef MLAPI_LC' that must be kept in place.
#define MLAPI_LC

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Teuchos includes
#ifdef HAVE_TEUCHOS
#include "Teuchos_PythonParameter.h"
#endif

// Epetra includes
#ifdef HAVE_EPETRA
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"

// Epetra python includes
#include "NumPyImporter.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#endif

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "MLAPI.h"
#include "MLAPI_PyMatrix.h"

%}

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include ML documentation
%include "ML_dox.i"

// General ignore directives
%ignore *::operator=;
%ignore *::operator[];

// External Trilinos package imports
#ifdef HAVE_TEUCHOS
%import "Teuchos.i"
#endif
#ifdef HAVE_EPETRA
%include "Epetra_RowMatrix_Utils.i"
%import  "Epetra.i"
#endif

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

///////////////////////
// ml_config support //
///////////////////////
%include "ml_config.h"

/////////////////////////////////////////
// ml_MultiLevelPreconditioner support //
/////////////////////////////////////////
%include "ml_MultiLevelPreconditioner.h"
#ifdef HAVE_TEUCHOS
#ifdef HAVE_EPETRA
namespace ML_Epetra
{
%extend MultiLevelPreconditioner
{
  // This is for MatrixPortal
  int SetParameterListAndNullSpace(PyObject* obj,
				   Epetra_MultiVector& NullSpace) {
    Teuchos::ParameterList * List = Teuchos::pyDictToNewParameterList(obj);
    if (List == NULL) List = new Teuchos::ParameterList();
    // WARNING: THIS IS DELICATE, NULLSPACE SHOULD NOT DISAPPEAR
    // otherwise the pointer here stored will vanish. This function should
    // be used only through the MatrixPortal
    double* NullSpacePtr = (double*)NullSpace.Values();
    int NullSpaceDim = NullSpace.NumVectors();
    List->set("null space: type", "pre-computed");
    List->set("null space: vectors", NullSpacePtr);
    List->set("null space: dimension", NullSpaceDim);
    self->SetParameterList(*List);
    delete List;
    return(0);
  }
}
}
#endif
#endif

/////////////////////////////
// MLAPI_Workspace support //
/////////////////////////////
%include "MLAPI_Workspace.h"

//////////////////////////////
// MLAPI_BaseObject support //
//////////////////////////////
%include "MLAPI_BaseObject.h"
namespace MLAPI
{
%extend BaseObject
{
  // Define the __str__() method, used by the python str() operator
  // on any object given to the python print command.
  std::string __str__()
  {
    std::stringstream os;
    self->Print(os);                  // Put the output in os
    std::string s = os.str();         // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}
}

//////////////////////////////
// MLAPI_CompObject support //
//////////////////////////////
%include "MLAPI_CompObject.h"

//////////////////////////////
// MLAPI_TimeObject support //
//////////////////////////////
%include "MLAPI_TimeObject.h"

////////////////////////////////
// MLAPI_BaseOperator support //
////////////////////////////////
%warnfilter(473)     MLAPI::BaseOperator;
%feature("director") MLAPI::BaseOperator;
%include "MLAPI_BaseOperator.h"

/////////////////////////
// MLAPI_Space support //
/////////////////////////
%include "MLAPI_Space.h"
namespace MLAPI
{
%extend Space
{
  PyObject*  GetMyGlobalElements()
  {
    int n = self->GetNumMyElements();
    if (n <= 0)
      goto fail;
    PyObject* MyGlobalElements_Python,* item;
    MyGlobalElements_Python = PyList_New(n);
    for (int i = 0 ; i < n ; ++i)
    {
      item = PyInt_FromLong((*self)(i));
      PyList_SetItem(MyGlobalElements_Python, i, item);
    }
    return(MyGlobalElements_Python);
  fail:
    Py_INCREF(Py_None);
    return Py_None;
  }
}
}

///////////////////////////////
// MLAPI_MultiVector support //
///////////////////////////////
%include "MLAPI_MultiVector.h"
namespace MLAPI
{
%extend MultiVector
{
  MultiVector __add__(MultiVector& rhs)
  {
    MultiVector res = *self + rhs;
    return(res);
  }
  MultiVector __sub__(MultiVector& rhs)
  {
    return(*self - rhs);
  }
  double __mul__(MultiVector& rhs)
  {
    return(*self * rhs);
  }
  MultiVector __mul__(double rhs)
  {
    return(*self * rhs);
  }
  PyObject * __setitem__(PyObject* args, double value)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row, Col;
    if (!(PyArg_ParseTuple(args, "OO:MLAPI_MultiVector___setitem__",
			   &rowObj, &colObj)      &&
	  SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	  SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))   ))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    (*self)(Row, Col) = value;
    return Py_BuildValue("");
  }
  PyObject * __getitem__(PyObject* args)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    if (!(PyArg_ParseTuple(args, "OO:MLAPI_MultiVector___getitem__",
			   &rowObj, &colObj)      &&
	  SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	  SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))    ))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    return SWIG_From_double((*self)(Row, Col));
  }
}
}

/////////////////////////////////////
// MLAPI_MultiVector_Utils support //
/////////////////////////////////////
%include "MLAPI_MultiVector_Utils.h"

////////////////////////////
// MLAPI_Operator support //
////////////////////////////
%include "MLAPI_Operator.h"
#ifdef HAVE_EPETRA
namespace MLAPI
{
%extend Operator
{
  PyObject* __getitem__(PyObject* args)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    Epetra_RowMatrix* Matrix = self->GetRCPRowMatrix().get();
    if (SWIG_IsOK(SWIG_AsVal_int(args, &Row)))
    {
      return Epetra_RowMatrix_GetEntries(*Matrix, Row);
    }
    else if (PyArg_ParseTuple(args, "OO:MLAPI_Operator___getitem__",
			      &rowObj, &colObj)      &&
	     SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	     SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))    )
    {
      return Epetra_RowMatrix_GetEntry(*Matrix, Row, Col);
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      return NULL;
    }
  }
  MultiVector __mul__(MultiVector& rhs)
  {
    return(*self * rhs);
  }
  Operator __add__(Operator& rhs)
  {
    return(*self + rhs);
  }
  Operator __sub__(Operator& rhs)
  {
    return(*self - rhs);
  }
  Operator __mul__(Operator& rhs)
  {
    if (self == &rhs)
      return(*self * Duplicate(rhs));
    else
      return(*self * rhs);
  }
  Operator __mul__(double rhs)
  {
    return(*self * rhs);
  }
  Operator __div__(double rhs)
  {
    return(*self * (1.0 / rhs));
  }
}
}
#endif

///////////////////////////////////
// MLAPI_InverseOperator support //
///////////////////////////////////
%include "MLAPI_InverseOperator.h"
#ifdef HAVE_TEUCHOS
namespace MLAPI
{
  %extend InverseOperator
  {
    MultiVector __mul__(MultiVector& rhs)
    {
      return(*self * rhs);
    }
    bool Reshape(const Operator& Op, const std::string Type, PyObject* obj)
    {
      Teuchos::ParameterList * List = Teuchos::pyDictToNewParameterList(obj);
      if (List == NULL) return(false);
      else
	self->Reshape(Op, Type, *List);
      delete List;
      return(true);
    }
  }
}
#endif

//////////////////////////////////
// MLAPI_Operator_Utils support //
//////////////////////////////////
%include "MLAPI_Operator_Utils.h"

//////////////////////////////////////
// MLAPI_EpetraBaseOperator support //
//////////////////////////////////////
#ifdef HAVE_EPETRA
%include "MLAPI_EpetraBaseOperator.h"
#endif

//////////////////////////
// MLAPI_Krylov support //
//////////////////////////
%include "MLAPI_Krylov.h"

///////////////////////////////
// MLAPI_Expressions support //
///////////////////////////////
%include "MLAPI_Expressions.h"

///////////////////////////////
// MLAPI_Aggregation support //
///////////////////////////////
%include "MLAPI_Aggregation.h"

///////////////////////
// MLAPI_Eig support //
///////////////////////
%include "MLAPI_Eig.h"

///////////////////////////
// MLAPI_Gallery support //
///////////////////////////
%include "MLAPI_Gallery.h"

////////////////////////////
// MLAPI_PyMatrix support //
////////////////////////////
#ifdef HAVE_EPETRA
%include "MLAPI_PyMatrix.h"
%extend PyMatrix
{
  PyObject * __setitem__(PyObject* args, double val)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    if (!(PyArg_ParseTuple(args, "OO:PyMatrix___setitem__",
			   &rowObj, &colObj)      &&
	  SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	  SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))   ))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    self->SetElement(Row, Col, val);
    return Py_BuildValue("");
  }
  PyObject* __getitem__(PyObject* args)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    if (SWIG_IsOK(SWIG_AsVal_int(args, &Row)))
    {
      return Epetra_RowMatrix_GetEntries(*(self->GetMatrix()), Row);
    }
    else if (PyArg_ParseTuple(args, "OO:PyMatrix___getitem__",
			      &rowObj, &colObj)      &&
	     SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	     SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))    )
    {
      return Epetra_RowMatrix_GetEntry(*(self->GetMatrix()), Row, Col);
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      return NULL;
    }
  }
}
#endif

#ifdef HAVE_TEUCHOS
%inline
{
  MLAPI::Operator GetPNonSmoothed(const MLAPI::Operator& A,
				  const MLAPI::MultiVector& ThisNS,
				  MLAPI::MultiVector& NextNS,
				  PyObject* obj)
  {
    Teuchos::ParameterList * List = Teuchos::pyDictToNewParameterList(obj);
    MLAPI::Operator Ptent;
    MLAPI::GetPtent(A, *List, ThisNS, Ptent, NextNS);
    delete List;
    return(Ptent);
  }

  bool Iterate(const MLAPI::Operator& A, const MLAPI::MultiVector& LHS,
	       const MLAPI::MultiVector& RHS, const MLAPI::BaseOperator& Prec, 
	       PyObject* obj)
  {
    Teuchos::ParameterList * List = Teuchos::pyDictToNewParameterList(obj);
    if (List == NULL) return(false);
    Krylov(A, LHS, RHS, Prec, *List);
    delete List;
    return(true);
  }
}
#endif

%pythoncode
%{
# Call MPI_Init if appropriate
Init()

# Arrange for MPI_Finalize to be called at exit, if appropriate
import atexit
atexit.register(Finalize)
%}

// Turn off the exception handling
%exception;
