// ///////////////////////////////////////////////////////////////////
// TSFCoreLinearOpTesterDecl.hpp

#ifndef TSFCORE_LINEAR_OP_TESTER_DECL_HPP
#define TSFCORE_LINEAR_OP_TESTER_DECL_HPP

#include "TSFCoreTypes.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace TSFCore {

///
/** Testing class for <tt>LinearOp</tt>.
 *
 */
template<class Scalar>
class LinearOpTester {
public:

  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, warning_tol );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, error_tol );

	///
	LinearOpTester(
		const double      warning_tol = 1e-13
		,const double     error_tol   = 1e-10
		);

	///
	bool check(
		const LinearOp<Scalar>      &op
		,std::ostream               *out
		) const;

}; // class LinearOpTester

} // namespace TSFCore

#endif // TSFCORE_LINEAR_OP_TESTER_DECL_HPP
