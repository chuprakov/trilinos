// ////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.h

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H

#include "MeritFunc_PenaltyParamUpdate_AddedStep.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for MeritFunc_PenaltyParamUpdate_AddedStep from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group MeritFuncPenaltyParamUpdate {
	    small_mu     = 1e-6;
	    min_mu_ratio = 1e-8;
	    mult_factor  = 1e-4;
	    kkt_near_sol = 1.0;
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[small_mu] The smallest mu allows when away from the soltion.\\
  *		Example: small_mu = 1e-6;
  *	\item[min_mu_ratio] This bounds the smallest mu(i) as:\\
  *		min(mu(i))/max(mu(i)) >= min_mu_ratio.\\
  *		Example: min_mu_ratio = 1e-4;
  *	\item[mult_factor] Multiplicative factor for mu(j) = (1.0+mult_factor) * abs(lambda(j)).\\
  *		Example: mult_factor = 1e-4;
  *	\item[kkt_near_sol] When the total kkt_error is below kkt_near_sol a safer
  *		penalty update will be used.\\
  *		Example: kkt_near_sol = 1.0;
  *	\end{description}
  */
class MeritFunc_PenaltyParamUpdate_AddedStepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			MeritFunc_PenaltyParamUpdate_AddedStep >
{
public:

	///
	MeritFunc_PenaltyParamUpdate_AddedStepSetOptions(
		MeritFunc_PenaltyParamUpdate_AddedStep* target = 0 );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class MeritFunc_PenaltyParamUpdate_AddedStepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H
