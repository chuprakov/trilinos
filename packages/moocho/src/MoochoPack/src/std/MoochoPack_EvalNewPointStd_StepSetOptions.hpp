// ////////////////////////////////////////////////////////////////
// MoochoPack_EvalNewPointStd_StepSetOptions.hpp
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

#ifndef EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H
#define EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H

#include "MoochoPack_EvalNewPointStd_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for EvalNewPointStd_Step from an \c OptionsFromStream object.
 *
 * The default options group name is EvalNewPointStd.
 *
 * The options group is:
 *
 \verbatim
    options_group EvalNewPointStd {
        fd_deriv_testing = FD_DEFAULT;
        decomp_sys_teting = DST_DEFAULT;
        decomp_sys_teting_print_level = DSPL_USE_GLOBAL;
    }
 \verbatim
 *
 * <ul>
 * <li> <b>fd_deriv_testing</b>: Determines if finite differerece testing of the 
 *      derivatives of the Gc and Gf.  See the class \c EvalNewPointStd_Step
 *      and its printed algorithm for more details.
 *      <ul>
 *      <li> <b>FD_DEFAULT</b>: The global flag check_results determines
 *           if the tests are performed.
 *      <li> <b>FD_TEST</b>: The tests are performed reguardless the
 *           value of check_results
 *      <li> <b>FD_NO_TEST</b>: The tests are not performed reguardless the
 *           value of check_results
 *      </ul>
 * <li> <b>decomp_sys_testing</b>: Determines if the range/null decomposition of
 *      Gc and Gh is performed.  See the class \c EvalNewPointStd_Step
 *      and its printed algorithm for more details.
 *      <ul>
 *      <li> <b>DST_DEFAULT</b>: The global flag check_results determines
 *           if the tests are performed.
 *      <li> <b>DST_TEST</b>: The tests are performed reguardless the
 *           value of check_results
 *      <li> <b>DST_NO_TEST</b>: The tests are not performed reguardless the
 *           value of check_results
 *      </ul>
 * </ul>
 */
class EvalNewPointStd_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			EvalNewPointStd_Step >
{
public:

	///
	EvalNewPointStd_StepSetOptions(
		  EvalNewPointStd_Step* target = 0
		, const char opt_grp_name[] = "EvalNewPointStd" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class EvalNewPointStd_StepSetOptions

}	// end namespace MoochoPack

#endif	// EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H
