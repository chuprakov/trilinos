// ////////////////////////////////////////////////////////////////
// QPSolverRelaxedTesterSetOptions.h

#ifndef QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H
#define QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H

#include "QPSolverRelaxedTester.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ConstrainedOptimizationPack {

///
/** Set options for QPSolverRelaxedTester from an
  * OptionsFromStream object.
  *
  * The default options group name is QPSolverRelaxedTester.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group QPSolverRelaxedTester {
	    opt_warning_tol   = 1e-10;  *** Tolerances for optimality conditions
	    opt_error_tol     = 1e-5;
	    feas_warning_tol  = 1e-10;  *** Tolerances for feasibility
	    feas_error_tol    = 1e-5;
	    comp_warning_tol  = 1e-10;  *** Tolerances for complementarity
	    comp_error_tol    = 1e-5;
	}
  \end{verbatim}
  */
class QPSolverRelaxedTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			QPSolverRelaxedTester >
{
public:

	///
	QPSolverRelaxedTesterSetOptions(
		  QPSolverRelaxedTester* target = 0
		, const char opt_grp_name[] = "QPSolverRelaxedTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class QPSolverRelaxedTesterSetOptions

}	// end namespace ConstrainedOptimizationPack

#endif	// QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H
