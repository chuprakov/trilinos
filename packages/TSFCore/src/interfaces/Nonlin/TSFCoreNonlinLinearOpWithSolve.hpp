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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearOpWithSolve.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_HPP

#include "TSFCoreNonlinLinearOpWithSolveDecl.hpp"
#include "TSFCoreNonlinLinearSolveOp.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVectorStdOps.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpWithSolve<Scalar> >
LinearOpWithSolve<Scalar>::clone_lows() const
{
	return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
LinearOpWithSolve<Scalar>::preconditioner() const
{
	return Teuchos::null;
}

// Overridden methods from LinearOp

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
LinearOpWithSolve<Scalar>::clone() const
{
	return this->clone_lows();
}

// Overridden methods from LinearSolveOp

template<class Scalar>
Teuchos::RefCountPtr<const LinearSolveOp<Scalar> >
LinearOpWithSolve<Scalar>::clone_lso() const
{
	return this->clone_lows();
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_HPP
