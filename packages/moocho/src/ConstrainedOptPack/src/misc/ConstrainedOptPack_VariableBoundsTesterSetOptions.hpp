// ////////////////////////////////////////////////////////////////
// VariableBoundsTesterSetOptions.h

#ifndef VARIABLE_BOUNDS_TESTER_SET_OPTIONS_H
#define VARIABLE_BOUNDS_TESTER_SET_OPTIONS_H

#include "VariableBoundsTester.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ConstrainedOptimizationPack {

///
/** Set options for \Ref{VariableBoundsTester} from an
  * OptionsFromStream object.
  *
  * The default options group name is VariableBoundsTester.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group VariableBoundsTester {
	    warning_tol   = 1e-10;
	    error_tol     = 1e-5;
	}
  \end{verbatim}
  */
class VariableBoundsTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			VariableBoundsTester >
{
public:

	///
	VariableBoundsTesterSetOptions(
		  VariableBoundsTester* target = 0
		, const char opt_grp_name[] = "VariableBoundsTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class VariableBoundsTesterSetOptions

}	// end namespace ConstrainedOptimizationPack

#endif	// VARIABLE_BOUNDS_TESTER_SET_OPTIONS_H
