// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_H

#include "ReducedHessianSecantUpdateBFGSProjected_Strategy.h"
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
	options_group ReducedHessianSecantUpdateBFGSProjected {
		proj_start_act_set_frac   = 0.8;   *** (+dbl)
        super_basic_mult_drop_tol = 1e-5;  *** (+dbl)
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[proj_start_act_set_frac] ToDo : Finish.
  * \item[super_basic_mult_drop_tol] ToDo: Finish.
  *	\end{description}
  */
class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			ReducedHessianSecantUpdateBFGSProjected_Strategy >
{
public:

	///
	ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions(
		ReducedHessianSecantUpdateBFGSProjected_Strategy* target = 0
		, const char opt_grp_name[] = "ReducedHessianSecantUpdateBFGSProjected" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_H
