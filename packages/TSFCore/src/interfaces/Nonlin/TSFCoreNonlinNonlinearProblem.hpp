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

// ////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblem.hpp

#ifndef TSFCORE_NONLINEAR_PROBLEM_HPP
#define TSFCORE_NONLINEAR_PROBLEM_HPP

#include "TSFCoreNonlinNonlinearProblemDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Nonlin {

//
// static members
//

template<class Scalar>
Scalar NonlinearProblem<Scalar>::infiniteBound()
{
	return 1e+50;
//	numeric_limits doesn't exist on GNU 2.96 compilers (HKT, 09/22/03). 
//	Returned large number from old version that was deemed large enough.
//	return std::numeric_limits<Scalar>::max();
//	To Do (HKT) : Create work-around for numeric_limits.
}

//
// object members
//

// Basic information

template<class Scalar>
int NonlinearProblem<Scalar>::Nu() const
{
	return 0;
}

template<class Scalar>
int NonlinearProblem<Scalar>::numResponseFunctions() const
{
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > space_g = this->space_g();
	return space_g.get() ? space_g->dim() : 0;
}

// VectorSpaces

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> > 
NonlinearProblem<Scalar>::space_u(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::space_u(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
NonlinearProblem<Scalar>::space_g() const
{
	return Teuchos::null;
}

// Bounds

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::uL(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::uL(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return uL(l);
}

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::uU(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::uU(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return uU(l);
}

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::gL() const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::gL(l): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
	return gL();
}

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::gU() const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::gL(l): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
	return gU();
}

// Initial values (guesses) for state and auxiliary variables

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::u0(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::uU(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return u0(l);
}

// Set and access calculation storage

template<class Scalar>
const Vector<Scalar>* NonlinearProblem<Scalar>::get_c() const
{
	return const_cast<NonlinearProblem*>(this)->get_c();
}

template<class Scalar>
void NonlinearProblem<Scalar>::set_g(Vector<Scalar>* g)
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::set_g(g): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
}

template<class Scalar>
Vector<Scalar>* NonlinearProblem<Scalar>::get_g()
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::get_g(): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
	return NULL;
}

template<class Scalar>
const Vector<Scalar>* NonlinearProblem<Scalar>::get_g() const
{
	return const_cast<NonlinearProblem*>(this)->get_g();
}

// Calculation methods

template<class Scalar>
void NonlinearProblem<Scalar>::calc_g(
	const Vector<Scalar>     &y
	,const Vector<Scalar>*   u[]
	,bool                    newPoint
	) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::calc_g(...): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
}

template<class Scalar>
void NonlinearProblem<Scalar>::reportFinalSolution(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    solved
		)
{
	// By default we just ignore the solution.
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLINEAR_PROBLEM_HPP
