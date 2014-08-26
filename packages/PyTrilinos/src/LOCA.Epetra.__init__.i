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

%define %loca_epetra_docstring
"
PyTrilinos.LOCA.Epetra is the python interface to namespace Epetra of
the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Epetra is to provide an extension of the
NOX.Epetra.Group to LOCA.  The python version of LOCA.Epetra supports
the following classes:

    * Group  - Extension of the NOX.Epetra.Group to LOCA
"
%enddef

%module(package      = "PyTrilinos.LOCA.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_epetra_docstring) __init__

%{
// System includes
#include <vector>

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

// PyTrilinos includes
#include "PyTrilinos_Teuchos_Util.h"
#include "PyTrilinos_Epetra_Util.h"

// Local Epetra includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"

// Epetra includes
#include "Epetra_DLLExportMacro.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Time.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

// NOX includes
#include "NOX.H"
#include "NOX_Epetra_Group.H"

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"
#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
#include "LOCA_Epetra.H"
#include "LOCA_Epetra_Group.H"

// Namespace flattening
using Teuchos::RCP;
using Teuchos::rcp;
%}

%ignore *::operator=;

// SWIG library includes
%include "stl.i"

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

%include "Epetra_DLLExportMacro.h"

// Teuchos and Epetra interface support
%import "Teuchos.i"
%include "Epetra_Base.i"    // For PyExc_EpetraError
%import "Epetra.i"

// Teuchos RCP support
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterface)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterfaceMVDX)
%teuchos_rcp(LOCA::PhaseTransition::AbstractGroup)
%teuchos_rcp(LOCA::TimeDependent::AbstractGroup)
%teuchos_rcp(LOCA::BorderedSystem::AbstractGroup)
%teuchos_rcp(LOCA::Homotopy::AbstractGroup)
%teuchos_rcp(LOCA::Homotopy::Group)
%teuchos_rcp(LOCA::Homotopy::DeflatedGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::ExtendedGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::Constraint)
%teuchos_rcp(LOCA::Pitchfork::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Pitchfork::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Abstract::Group)
%teuchos_rcp(LOCA::Abstract::TransposeSolveGroup)

// NOX interface support
%import "NOX.Abstract.i"
%import "NOX.Epetra.__init__.i"
%import "NOX.Epetra.Interface.i"

// Allow import from the parent directory
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
from .. import Abstract
%}

// LOCA base classes
%import(module="Extended") "LOCA_Extended_MultiAbstractGroup.H"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_ConstraintInterface.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
%import(module="PhaseTransition") "LOCA_PhaseTransition_AbstractGroup.H"
%import(module="TimeDependent") "LOCA_TimeDependent_AbstractGroup.H"
%import(module="BorderedSystem") "LOCA_BorderedSystem_AbstractGroup.H"
%import(module="Homotopy") "LOCA_Homotopy_AbstractGroup.H"
%import(module="Homotopy") "LOCA_Homotopy_Group.H"
%import(module="Homotopy") "LOCA_Homotopy_DeflatedGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%import(module="TurningPoint.MinimallyAugmented") "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
%import(module="TurningPoint.MinimallyAugmented") "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
%import(module="Hopf.MooreSpence") "LOCA_Hopf_MooreSpence_AbstractGroup.H"
%import(module="Hopf.MooreSpence") "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_Constraint.H"
%import(module="Pitchfork.MooreSpence") "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
%import(module="Pitchfork.MinimallyAugmented") "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
%import(module="Abstract") "LOCA_Abstract_Group.H"
%import(module="Abstract") "LOCA_Abstract_TransposeSolveGroup.H"

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
  catch(int errCode)
  {
    PyErr_Format(PyExc_EpetraError, "Error code = %d\nSee stderr for details", errCode);
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

/////////////////////////
// LOCA Epetra support //
/////////////////////////

#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
%include "LOCA_Epetra.H"

//////////////////////////////
// LOCA.Epetra.Group support //
//////////////////////////////

// temporarily ignore conflict-causing constructor.  TODO: fix this issue
%ignore LOCA::Epetra::Group::Group(Teuchos::RCP< LOCA::GlobalData > const &,Teuchos::ParameterList &,Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree > const &,NOX::Epetra::Vector &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,LOCA::ParameterVector const &);

%teuchos_rcp(LOCA::Epetra::Group)
%include "LOCA_Epetra_Group.H"
