// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearOpWithSolve.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_HPP

#include <assert.h>

#include "TSFCoreNonlinLinearOpWithSolveDecl.hpp"
#include "TSFCoreNonlinLinearSolveOp.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVectorStdOps.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearOpWithSolve<Scalar> >
LinearOpWithSolve<Scalar>::clone_lows() const
{
	return MemMngPack::null;
}

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearOp<Scalar> >
LinearOpWithSolve<Scalar>::preconditioner() const
{
	return MemMngPack::null;
}

// Overridden methods from LinearOp

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearOp<Scalar> >
LinearOpWithSolve<Scalar>::clone() const
{
	return this->clone_lows();
}

// Overridden methods from LinearSolveOp

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearSolveOp<Scalar> >
LinearOpWithSolve<Scalar>::clone_lso() const
{
	return this->clone_lows();
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_HPP
