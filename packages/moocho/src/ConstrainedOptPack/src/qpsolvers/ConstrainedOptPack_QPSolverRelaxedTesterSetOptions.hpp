// ////////////////////////////////////////////////////////////////
// ConstrainedOptPack_QPSolverRelaxedTesterSetOptions.hpp
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

#ifndef QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H
#define QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H

#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

///
/** Set options for QPSolverRelaxedTester from an
 * OptionsFromStream object.
 *
 * The default options group name is QPSolverRelaxedTester.
 *
 * The options group is:
 \verbatim
	options_group QPSolverRelaxedTester {
	    opt_warning_tol   = 1e-10;  *** Tolerances for optimality conditions
	    opt_error_tol     = 1e-5;
	    feas_warning_tol  = 1e-10;  *** Tolerances for feasibility
	    feas_error_tol    = 1e-5;
	    comp_warning_tol  = 1e-10;  *** Tolerances for complementarity
	    comp_error_tol    = 1e-5;
	}
  \endverbatim
  */
class QPSolverRelaxedTesterSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			QPSolverRelaxedTester >
{
public:

	///
	QPSolverRelaxedTesterSetOptions(
		  QPSolverRelaxedTester* target = 0
		, const char opt_grp_name[] = "QPSolverRelaxedTester" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class QPSolverRelaxedTesterSetOptions

}	// end namespace ConstrainedOptPack

#endif	// QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H
