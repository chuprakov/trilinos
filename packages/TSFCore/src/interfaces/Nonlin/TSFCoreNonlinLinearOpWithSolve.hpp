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
