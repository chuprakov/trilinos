// ////////////////////////////////////////////////////////////////
// LineSearchWatchDog_StepSetOptions.h

#ifndef LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
#define LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H

#include "LineSearchWatchDog_Step.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for \Ref{LineSearchWatchDog_Step} from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group LineSearchWatchDog {
		opt_kkt_err_threshold	= 1e-3; *** (+dbl)
		feas_kkt_err_threshold	= 1e-3; *** (+dbl)
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[opt_kkt_err_threshold] ToDo : Finish.
  *		Example: opt_kkt_err_threshold = 1e-1;
  *	\item[feas_kkt_err_threshold] ToDo : Finish.
  *		Example: feas_kkt_err_threshold = 1e-2;
  *	\end{description}
  */
class LineSearchWatchDog_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			LineSearchWatchDog_Step >
{
public:

	///
	LineSearchWatchDog_StepSetOptions(
		LineSearchWatchDog_Step* target = 0 );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class LineSearchWatchDog_StepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
