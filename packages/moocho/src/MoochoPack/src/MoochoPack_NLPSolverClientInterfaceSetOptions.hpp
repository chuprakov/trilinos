// ////////////////////////////////////////////////////////////////
// rSQPSolverClientInterfaceSetOptions.h

#ifndef RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H
#define RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H

#include "rSQPSolverClientInterface.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for rSQPSolverClientInterface from an
  * OptionsFromStream object.
  *
  * The default options group name is rSQPSolverClientInterface.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group rSQPSolverClientInterface {
        max_iter = ?;
        max_run_time = ?;  *** In minutes
        opt_tol = ?;
        feas_tol = ?;
        step_tol = ?;
        max_var_bounds_viol = ?;
		journal_output_level = PRINT_NOTHING;
	}
  \end{verbatim}
  *
  * See the class \Ref{rSQPSolverClientInterface} for a description of these
  * parameters.
  */
class rSQPSolverClientInterfaceSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			rSQPSolverClientInterface >
{
public:

	///
	rSQPSolverClientInterfaceSetOptions(
		  rSQPSolverClientInterface* target = 0
		, const char opt_grp_name[] = "rSQPSolverClientInterface" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class rSQPSolverClientInterfaceSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H