// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversIterativeLinearSolver.hpp

#ifndef TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP
#define TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP

#include "TSFCoreSolversIterativeLinearSolverDecl.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
MemMngPack::ref_count_ptr<const IterativeLinearSolver<Scalar> >
IterativeLinearSolver<Scalar>::clone() const
{
	return MemMngPack::null;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_ITERATIVE_LINEAR_SOLVER_HPP
