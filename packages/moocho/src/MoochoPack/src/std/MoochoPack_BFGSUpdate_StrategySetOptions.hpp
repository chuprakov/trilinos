// ////////////////////////////////////////////////////////////////
// MoochoPack_BFGSUpdate_StrategySetOptions.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef REDUCED_HESSIAN_BFGS_STD_STEP_SET_OPTIONS_H
#define REDUCED_HESSIAN_BFGS_STD_STEP_SET_OPTIONS_H

#include "MoochoPack_BFGSUpdate_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for BFGSUpdate_Strategy from an OptionsFromStream
  * object.
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
	void setOption( int option_num, const std::string& option_value );

};	// end class BFGSUpdate_StrategySetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_BFGS_STD_STEP_SET_OPTIONS_H
