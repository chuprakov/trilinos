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
