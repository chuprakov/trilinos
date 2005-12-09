// ////////////////////////////////////////////////////////////////
// MoochoPack_ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.hpp
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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_ReducedHessianSecantUpdateBFGSProjected_Strategy.hpp"
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
	options_group ReducedHessianSecantUpdateBFGSProjected {
		act_set_frac_proj_start   = 0.8;   *** (+dbl)
		project_error_tol         = 1e-5;  *** (+dbl) [0.0, 1.0]
        super_basic_mult_drop_tol = 1e-5;  *** (+dbl)
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[act_set_frac_proj_start] ToDo : Finish.
  *	\item[project_error_tol] ToDo : Finish.
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
		, const char opt_grp_name[] = "ReducedHessianSecantUpdatePBFGS" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H
