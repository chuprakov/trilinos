// ////////////////////////////////////////////////////////////////
// BFGSUpdate_StrategySetOptions.h

#ifndef REDUCED_HESSIAN_BFGS_STD_STEP_SET_OPTIONS_H
#define REDUCED_HESSIAN_BFGS_STD_STEP_SET_OPTIONS_H

#include "BFGSUpdate_Strategy.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for \Ref{BFGSUpdate_Strategy} from an
  * OptionsFromStream object.
  *
  * The default options group name is BFGSUpdate.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group BFGSUpdate {
	    rescale_init_identity   = true;
	    use_dampening           = true;
	    secant_testing          = DEFAULT;
	*    secant_testing          = TEST;
	*    secant_testing          = NO_TEST;
	    secant_warning_tol      = 1e-6;
	    secant_error_tol        = 1e-2;
	}
  \end{verbatim}
  */
class BFGSUpdate_StrategySetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			BFGSUpdate_Strategy >
{
public:

	///
	BFGSUpdate_StrategySetOptions(
		  BFGSUpdate_Strategy* target = 0
		, const char opt_grp_name[] = "BFGSUpdate" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class BFGSUpdate_StrategySetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// REDUCED_HESSIAN_BFGS_STD_STEP_SET_OPTIONS_H
