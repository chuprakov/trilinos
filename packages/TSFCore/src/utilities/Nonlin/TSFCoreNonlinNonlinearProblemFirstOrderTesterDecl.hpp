// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemFirstOrderTesterDecl.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_DECL_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_DECL_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreNonlinTypes.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Unit testing class for <tt>NonlinearProblemFirstOrder</tt>.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NonlinearProblemFirstOrderTester {
public:

	///
	enum EOutputLevel { BASIC_OUTPUT, VERBOSE_OUTPUT };

	///
	/** Test a <tt>NonlinearProblemFirstOrder</tt> object.
	 *
	 * ToDo: Finish documentation!
	 */
	bool doTest( NonlinearProblemFirstOrder<Scalar>* prob, std::ostream* out = NULL, EOutputLevel outputLevel = BASIC_OUTPUT ) const;

}; // NonlinearProblemFirstOrderTester

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_TESTER_DECL_HPP
