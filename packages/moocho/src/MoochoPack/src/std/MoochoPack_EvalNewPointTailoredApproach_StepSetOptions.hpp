// ////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproach_StepSetOptions.h

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H

#include "EvalNewPointTailoredApproach_Step.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for EvalNewPointTailoredApproach_Step from an
  * OptionsFromStream object.
  *
  * The default options group name is EvalNewPointTailoredApproach.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group EvalNewPointTailoredApproach {
		fd_deriv_testing   = FD_DEFAULT;
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[fd_deriv_testing] Determines if finite differerece testing of the 
  *		derivatives of the Gc and Gf.  See the class \Ref{EvalNewPointTailoredApproach_Step}
  *		and its printed algorithm for more details).
  *		\begin{description}
  *		\item[FD_DEFAULT]			The global flag check_results determines
  *									if the tests are performed.
  *		\item[FD_TEST]				The tests are performed reguardless the
  *									value of check_results
  *		\item[FD_NO_TEST]			The tests are not performed reguardless the
  *									value of check_results
  *		\end{description}
  */
class EvalNewPointTailoredApproach_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			EvalNewPointTailoredApproach_Step >
{
public:

	///
	EvalNewPointTailoredApproach_StepSetOptions(
		  EvalNewPointTailoredApproach_Step* target = 0
		, const char opt_grp_name[] = "EvalNewPointTailoredApproach" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class EvalNewPointTailoredApproach_StepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H
