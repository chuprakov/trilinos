// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H
#define REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H

#include "ReducedHessianSecantUpdateLPBFGS_Strategy.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for \Ref{ReducedHessianSecantUpdateBFGSProjected_Strategy} from a
  * OptionsFromStream object.
  *
  * The options group is (with the default name):
  *
  \begin{verbatim}
	options_group ReducedHessianSecantUpdateLPBFGS {
		min_num_updates_proj_start   = 0;      *** (+int)
		max_num_updates_proj_start   = 999999; *** (+int)
		num_superbasics_switch_dense = 500;    *** (+int)
		num_add_recent_updates       = 10;     *** (+int)
    }
  \end{verbatim}
  *
  * \begin{description}
  *	\item ToDo : Finish
  *	\end{description}
  */
class ReducedHessianSecantUpdateLPBFGS_StrategySetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
	, public OptionsFromStreamPack::SetOptionsToTargetBase<
	      ReducedHessianSecantUpdateLPBFGS_Strategy >
	{
public:

	///
	ReducedHessianSecantUpdateLPBFGS_StrategySetOptions(
		ReducedHessianSecantUpdateLPBFGS_Strategy* target = 0
		, const char opt_grp_name[] = "ReducedHessianSecantUpdateLPBFGS" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateLPBFGS_StrategySetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H
