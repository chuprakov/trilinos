// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemTesterDecl.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_TESTER_DECL_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_TESTER_DECL_HPP

#include <iosfwd>

#include "TSFCoreNonlinTypes.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Unit testing class for <tt>NonlinearProblem</tt>.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NonlinearProblemTester {
public:

	///
	enum EOutputLevel { BASIC_OUTPUT, VERBOSE_OUTPUT };

	///
	/** Test a <tt>NonlinearProblem</tt> object.
	 *
	 * ToDo: Finish documentation!
	 */
	bool doTest( NonlinearProblem<Scalar>* prob, std::ostream* out = NULL, EOutputLevel outputLevel = BASIC_OUTPUT ) const;

}; // NonlinearProblemTester

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_TESTER_DECL_HPP
