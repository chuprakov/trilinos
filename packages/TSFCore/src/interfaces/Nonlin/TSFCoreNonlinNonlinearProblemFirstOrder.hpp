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

// ////////////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemFirstOrder.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrderDecl.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
bool NonlinearProblemFirstOrder<Scalar>::adjointSupported() const
{
	return true;
}

// Factories for linear operators

template<class Scalar>
Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOp<Scalar> > >
NonlinearProblemFirstOrder<Scalar>::factory_DcDu(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::factory_DcDu(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return Teuchos::null;
}

// Transpose arguments

template<class Scalar>
ETransp NonlinearProblemFirstOrder<Scalar>::opDcDu(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::opDcDu(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return NOTRANS;
}

// Set and access calculation storage

template<class Scalar>
const LinearOpWithSolve<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DcDy() const
{
	return const_cast<NonlinearProblemFirstOrder*>(this)->get_DcDy();
}

template<class Scalar>
void NonlinearProblemFirstOrder<Scalar>::set_DcDu(int l, LinearOp<Scalar>* DcDu_l)
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::set_DcDu(l,DcDu_l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
}

template<class Scalar>
LinearOp<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DcDu(int l)
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::get_DcDu(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return NULL;
}

template<class Scalar>
const LinearOp<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DcDu(int l) const
{
	return const_cast<NonlinearProblemFirstOrder*>(this)->get_DcDu(l);
}

template<class Scalar>
void NonlinearProblemFirstOrder<Scalar>::set_DgDy(MultiVector<Scalar>* DgDy)
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::set_DgDy(DgDy): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
}

template<class Scalar>
MultiVector<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DgDy()
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::set_DgDy(): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
	return NULL;
}

template<class Scalar>
const MultiVector<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DgDy() const
{
	return const_cast<NonlinearProblemFirstOrder*>(this)->get_DgDy();
}

template<class Scalar>
void NonlinearProblemFirstOrder<Scalar>::set_DgDu(int l, MultiVector<Scalar>* DgDu_l)
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::set_DgDu(l,DgDu_l): Error, Must be overridden "
		"in subclass along with Nu() > 0 and space_g().get() != NULL!"
		);
}

template<class Scalar>
MultiVector<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DgDu(int l)
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::set_DgDu(l,DgDu_l): Error, Must be overridden "
		"in subclass along with Nu() > 0 and space_g().get() != NULL!"
		);
	return NULL;
}

template<class Scalar>
const MultiVector<Scalar>* NonlinearProblemFirstOrder<Scalar>::get_DgDu(int l) const
{
	return const_cast<NonlinearProblemFirstOrder*>(this)->get_DgDu(l);
}

// Calculation methods

template<class Scalar>
void NonlinearProblemFirstOrder<Scalar>::calc_DcDu(
	int                      l
	,const Vector<Scalar>    &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::clac_DcDu(...): Error, Must be overridden "
		"in subclass along with Nu() > 0!"
		);
}

template<class Scalar>
void NonlinearProblemFirstOrder<Scalar>::calc_DgDy(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::clac_DgDy(...): Error, Must be overridden "
		"in subclass along with space_g().get() != NULL!"
		);
}

template<class Scalar>
void NonlinearProblemFirstOrder<Scalar>::calc_DgDu(
	int                      l
	,const Vector<Scalar>    &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblemFirstOrder<Scalar>::clac_DgDy(...): Error, Must be overridden "
		"in subclass along with Nu() > 0 and space_g().get() != NULL!"
		);
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_HPP
