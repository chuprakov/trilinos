// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblemFirstOrder.cpp

#include "Epetra_NonlinearProblemFirstOrder.hpp"

#define EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED( FNC_NAME ) \
  TEST_FOR_EXCEPTION( \
    true, std::logic_error \
    ,"Error, the function " << (FNC_NAME) << " given a default implementation in " \
    "Epetra::NonlinearProblemFirstOrder is not overridden in the class " \
    << typeid(*this).name() << " or should not be called by the client!" \
  )

namespace Epetra {

// Adjoints supported?

bool NonlinearProblemFirstOrder::adjointSupported() const
{
  return true;
}

// Factories for linear operators

bool NonlinearProblemFirstOrder::use_DcDu_op(int l) const
{
  EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED("use_DcDu_op(l)");
  return false;
}

Teuchos::RefCountPtr<Epetra_Operator>
NonlinearProblemFirstOrder::create_DcDu_op(int l) const
{
  EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED("create_DcDu_op(l)");
  return Teuchos::null;
}

bool NonlinearProblemFirstOrder::DcDy_op_is_const() const
{
	return false;
}

bool NonlinearProblemFirstOrder::specialized_DcDy_prec() const
{
  return false;
}

Teuchos::RefCountPtr<Epetra_Operator>
NonlinearProblemFirstOrder::create_DcDy_prec() const
{
  return Teuchos::null;
}

bool NonlinearProblemFirstOrder::DcDy_prec_is_const() const
{
	return false;
}

Teuchos::RefCountPtr<Epetra_MultiVector>
NonlinearProblemFirstOrder::create_DcDu_mv(int l) const
{
  if(this->use_DcDu_op(l)) {
    return Teuchos::null;
  }
  else {
    return Teuchos::rcp(
      new Epetra_MultiVector(
        *this->map_c()
        ,map_u(l)->NumGlobalElements()
        )
      );
  }
}

bool NonlinearProblemFirstOrder::DcDu_is_const(int l) const
{
	return false;
}

// Transpose arguments

ETransp NonlinearProblemFirstOrder::opDcDu(int l) const
{
  TEST_FOR_EXCEPTION(
    this->use_DcDu_op(l), std::logic_error
    ,"NonlinearProblemFirstOrder::opDcDu(l): Error, the default implementation for "
    "this function is for the creation of Epetra_MultiVector objects for DcDu.  For "
    "different behavior, this function must be overridden by the subclass "
    << typeid(*this).name() << " or one of its base classes!" \
    );
  return NOTRANS;
}

// Calculation methods

void NonlinearProblemFirstOrder::calc_Dg(
  const Epetra_Vector           &y
  ,const Epetra_Vector*         u[]
  ,Epetra_Vector                *g
  ,Epetra_MultiVector           *DgDy
  ,Epetra_MultiVector*          DgDu[]
  ) const
{
  EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED("calc_Dg(...)");
}

// Overridden from NonlinearProblem

void NonlinearProblemFirstOrder::calc_c(
  const Epetra_Vector     &y
  ,const Epetra_Vector*   u[]
  ,Epetra_Vector          *c
  ) const
{
  calc_Dc( y, u, c, NULL, NULL, NULL );
}

void NonlinearProblemFirstOrder::calc_g(
  const Epetra_Vector     &y
  ,const Epetra_Vector*   u[]
  ,Epetra_Vector          *g
  ) const
{
  calc_Dg( y, u, g, NULL, NULL );
}

} // namespace Epetra
