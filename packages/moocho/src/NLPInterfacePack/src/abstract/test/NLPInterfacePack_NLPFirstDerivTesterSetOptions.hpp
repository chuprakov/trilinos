// ////////////////////////////////////////////////////////////////
// NLPFirstDerivativesTesterSetOptions.h

#ifndef NLP_FIRST_DERIVATIVES_TESTER_SET_OPTIONS_H
#define NLP_FIRST_DERIVATIVES_TESTER_SET_OPTIONS_H

#include "NLPFirstDerivativesTester.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace NLPInterfacePack {
namespace TestingPack {

///
/** Set options for NLPFirstDerivativesTester from an
  * OptionsFromStream object.
  *
  * The default options group name is NLPFirstDerivativesTester.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group NLPFirstDerivativesTester {
	    warning_tol   = 1e-6;
	    error_tol     = 1e-3;
	}
  \end{verbatim}
  */
class NLPFirstDerivativesTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPFirstDerivativesTester >
{
public:

	///
	NLPFirstDerivativesTesterSetOptions(
		  NLPFirstDerivativesTester* target = 0
		, const char opt_grp_name[] = "NLPFirstDerivativesTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class NLPFirstDerivativesTesterSetOptions

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_DERIVATIVES_TESTER_SET_OPTIONS_H