// ////////////////////////////////////////////////////////////////
// MoochoPack_ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.hpp
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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H
#define REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_ReducedHessianSecantUpdateLPBFGS_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for ReducedHessianSecantUpdateBFGSProjected_Strategy
  * from a OptionsFromStream object.
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
	void setOption( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateLPBFGS_StrategySetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H
