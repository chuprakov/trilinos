// ////////////////////////////////////////////////////////////////////////
// TSFCoreSolversTypes.hpp

#ifndef TSF_CORE_SOLVERS_TYPES_HPP
#define TSF_CORE_SOLVERS_TYPES_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {
namespace Solvers {

///
const int DEFAULT_MAX_ITER = -1;

///
enum ESolveReturn {
	SOLVED_TO_TOL          ///< The linear system(s) solved to tolerance
	,MAX_ITER_EXCEEDED     ///< <tt>max_iter</tt> exceeded
};

///
struct SolveReturn {
	SolveReturn(ESolveReturn solve_status_in, int num_iter_in) :solve_status(solve_status_in),num_iter(num_iter_in) {} 
	ESolveReturn    solve_status;
	int             num_iter;
};

namespace Exceptions {

/** \defgroup TSFCoreSolversExceptions_grp Basic TSFCore::Solvers exception types.
 */
//@{

/// Thrown if a linear solve failed
class FailureToConverge : public std::logic_error
{public: FailureToConverge(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a the operator is defective in some way (related to linear solver method)
class SolverBreakdown : public FailureToConverge
{public: SolverBreakdown(const std::string& what_arg) : FailureToConverge(what_arg) {}};

/// Thrown if operator turns out to be indefinite
class Indefinite : public SolverBreakdown
{public: Indefinite(const std::string& what_arg) : SolverBreakdown(what_arg) {}};

//@}

} // namespace Exceptions

template<class Scalar> class Norm;
template<class Scalar> class SolverState;
template<class Scalar> class ConvergenceTester;
template<class Scalar> class IterativeLinearSolver;

} // namespace Solvers
} // namespace TSFCore

#endif // TSF_CORE_SOLVERS_TYPES_HPP
