// ////////////////////////////////////////////////////////////////
// QPSolverRelaxedQPSchurSetOptions.h

#ifndef QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H
#define QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H

#include "QPSolverRelaxedQPSchur.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ConstrainedOptimizationPack {

///
/** Set options for QPSolverRelaxedQPSchur from an
  * OptionsFromStream object.
  *
  * The default options group name is QPSolverRelaxedQPSchur.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group QPSolverRelaxedQPSchur {
	*    max_qp_iter_frac  = 10.0;   *** (+dbl) max_qp_itr = max_qp_itr_frac * (# variables)
	*    inequality_pick_policy = ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY;
	*    inequality_pick_policy = ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY; *** not supported yet!
	*    inequality_pick_policy = ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY;
	*    bounds_tol        = 1e-10;  *** (+dbl) feasibility tolerance for bound constriants
	*    inequality_tol    = 1e-10;  *** (+dbl) feasibility tolerance for general inequality constriants
	*    equality_tol      = 1e-10;  *** (+dbl) feasibility tolerance for general equality constriants
	*    loose_feas_tol    = 1e-9;   *** (+dbl) (Expert use only)
	*    dual_infeas_tol   = 1e-12;  *** (+dbl) allowable dual infeasiblity before error
	*    huge_primal_step  = 1e+20;  *** (+dbl) value of a near infinite primal step
	*    huge_dual_step    = 1e+20;  *** (+dbl) value of a near infinite dual step
	*    bigM              = 1e+10;  *** (+dbl) value or relaxation penalty in objective
	*    warning_tol   = 1e-10;  *** Testing warning tolerance
	*    error_tol     = 1e-5;   *** Testing error tolerance
	*    print_level = USE_INPUT_ARG;  *** Use the input argument to solve_qp(...)
	*    print_level = NO_OUTPUT;
	*    print_level = OUTPUT_BASIC_INFO;
	*    print_level = OUTPUT_ITER_SUMMARY;
	*    print_level = OUTPUT_ITER_STEPS;
	*    print_level = OUTPUT_ACT_SET;
	*    print_level = OUTPUT_ITER_QUANTITIES;
	}
  \end{verbatim}
  */
class QPSolverRelaxedQPSchurSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			QPSolverRelaxedQPSchur >
{
public:

	///
	QPSolverRelaxedQPSchurSetOptions(
		  QPSolverRelaxedQPSchur* target = 0
		, const char opt_grp_name[] = "QPSolverRelaxedQPSchur" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class QPSolverRelaxedQPSchurSetOptions

}	// end namespace ConstrainedOptimizationPack

#endif	// QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H
