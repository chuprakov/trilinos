// ////////////////////////////////////////////////////////////////
// MoochoPack_NLPSolverClientInterfaceSetOptions.hpp
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

#ifndef RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H
#define RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H

#include "MoochoPack_NLPSolverClientInterface.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Set options for NLPSolverClientInterface from an \c OptionsFromStream object.
 *
 * The default options group name is NLPSolverClientInterface.
 *
 * The options group is:
 *
 \verbatim
	options_group NLPSolverClientInterface {
        max_iter = ?;
        max_run_time = ?;  *** In minutes
        opt_tol = ?;
        feas_tol = ?;
        step_tol = ?;
		journal_output_level = ?;
		journal_print_digits = ?;
		check_results = ?;
		calc_conditioning = ?
	}
 \endverbatim
 *
 * See the class \c NLPSolverClientInterface for a description of these
 * parameters.
 */
class NLPSolverClientInterfaceSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPSolverClientInterface >
{
public:

	///
	NLPSolverClientInterfaceSetOptions(
		  NLPSolverClientInterface* target = 0
		, const char opt_grp_name[] = "NLPSolverClientInterface" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class NLPSolverClientInterfaceSetOptions

}	// end namespace MoochoPack

#endif	// RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H
