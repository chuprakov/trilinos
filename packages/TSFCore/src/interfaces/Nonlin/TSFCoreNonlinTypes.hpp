// ///////////////////////////////////////////////////////////////
// TSFCoreNonlinTypes.hpp

#ifndef TSF_CORE_NONLIN_TYPES_HPP
#define TSF_CORE_NONLIN_TYPES_HPP

#include "TSFCoreSolversTypes.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Enumeration for determining if an operator is nonsingular or not.
 */
enum ENonsingStatus {
	OP_NONSINGULAR            ///< The operator is known to be (or should be) nonsingular
	,OP_SINGULAR              ///< The oeprator is known to be (or should be) singular
	,OP_SINGULARITY_UNKNOWN   ///< Unknown if the operator is singular or nonsingular
};

template<class Scalar> class LinearSolveOp;
template<class Scalar> class LinearOpWithSolve;
template<class Scalar> class NonlinearProblem;
template<class Scalar> class NonlinearProblemFirstOrder;

} // namespace Nonlin
} // namespace TSFCore

#endif // TSF_CORE_NONLIN_TYPES_HPP
