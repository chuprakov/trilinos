// ////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblem.hpp

#ifndef TSFCORE_NONLINEAR_PROBLEM_HPP
#define TSFCORE_NONLINEAR_PROBLEM_HPP

#include "TSFCoreNonlinNonlinearProblemDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "ThrowException.hpp"

namespace TSFCore {
namespace Nonlin {

//
// static members
//

template<class Scalar>
Scalar NonlinearProblem<Scalar>::infiniteBound()
{
	return std::numeric_limits<Scalar>::max();
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
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > space_g = this->space_g();
	return space_g.get() ? space_g->dim() : 0;
}

// VectorSpaces

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > 
NonlinearProblem<Scalar>::space_u(int l) const
{
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::space_u(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return MemMngPack::null;
}

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >
NonlinearProblem<Scalar>::space_g() const
{
	return MemMngPack::null;
}

// Bounds

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::uL(int l) const
{
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::uL(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return uL(l);
}

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::uU(int l) const
{
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::uU(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	return uU(l);
}

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::gL() const
{
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::gL(l): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
	return gL();
}

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::gU() const
{
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::gL(l): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
	return gU();
}

// Initial values (guesses) for state and auxiliary variables

template<class Scalar>
const Vector<Scalar>& NonlinearProblem<Scalar>::u0(int l) const
{
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::uU(l): Error, Must be overridden in subclass along with Nu() > 0!"
		);
	u0(l);
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
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::set_g(g): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
}

template<class Scalar>
Vector<Scalar>* NonlinearProblem<Scalar>::get_g()
{
	THROW_EXCEPTION(
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
	THROW_EXCEPTION(
		true,std::logic_error
		,"NonlinearProblem<Scalar>::calc_g(...): Error, Must be overridden in subclass along with space_g().get() != NULL!"
		);
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLINEAR_PROBLEM_HPP
