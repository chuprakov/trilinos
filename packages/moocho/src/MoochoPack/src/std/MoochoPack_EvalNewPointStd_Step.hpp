// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointStd_Step.hpp
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

#ifndef EVAL_NEW_POINT_STD_STEP_H
#define EVAL_NEW_POINT_STD_STEP_H

#include "DecompositionSystemHandler_Strategy.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "ConstrainedOptPack/src/decompositions/DecompositionSystemTester.hpp"
#include "ConstrainedOptPack/src/misc/VariableBoundsTester.hpp"
#include "NLPInterfacePack/src/abstract/test/NLPFirstDerivTester.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Standard new point evaluation step class.
 *
 * This class calculates \c Gc, \c Gh, updates the range/null decompositon matrices
 * \c Z, \c Y, \c R, \c Uz, \c Uy \c Vz and \c Vy and calculates \c Gf, \c c,
 * \c h, and \c f in that order.
 */
class EvalNewPointStd_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/** @name Public types */
	//@{

	///
	enum EFDDerivTesting   { FD_DEFAULT,  FD_TEST,  FD_NO_TEST  };

	//@}

	/** @name Constructors / initializers */
	//@{

	/// «std comp» members for range/null decomposition handler
	STANDARD_COMPOSITION_MEMBERS( DecompositionSystemHandler_Strategy, decomp_sys_handler )
	/// «std comp» members for first derivative tester object
	STANDARD_COMPOSITION_MEMBERS( NLPFirstDerivTester, deriv_tester )
	/// «std comp» Members for variable bounds tester object
	STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester )
	/// «std comp» members for decomp_sys tester tester object
	STANDARD_COMPOSITION_MEMBERS( DecompositionSystemTester, decomp_sys_tester )
	///
	/** Set how and if finite derivatives are tested.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDDerivTesting, fd_deriv_testing )
	///
	/** Set how and if the decomposition system is tested.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( DecompositionSystemHandler_Strategy::EDecompSysTesting, decomp_sys_testing )
	///
	/** Set how to set the print level for decomp_sys_tester (only if testing).
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( DecompositionSystemHandler_Strategy::EDecompSysPrintLevel, decomp_sys_testing_print_level )

	///
	/** Constructor.
	 *
	 * new_point == true by default.
	 */
	EvalNewPointStd_Step(
		const decomp_sys_handler_ptr_t                              &decomp_sys_handler
		,const deriv_tester_ptr_t                                   &deriv_tester
		,const bounds_tester_ptr_t                                  &bounds_tester
		,const decomp_sys_tester_ptr_t                              &decomp_sys_tester
		,EFDDerivTesting                                            fd_deriv_testing               = FD_DEFAULT
		,DecompositionSystemHandler_Strategy::EDecompSysTesting     decomp_sys_testing             = DecompositionSystemHandler_Strategy::DST_DEFAULT
		,DecompositionSystemHandler_Strategy::EDecompSysPrintLevel  decomp_sys_testing_print_level = DecompositionSystemHandler_Strategy::DSPL_USE_GLOBAL
		);

	//@}

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

private:

	// Not defined and not to be called
	EvalNewPointStd_Step();

};	// end class EvalNewPointStd_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_STD_STEP_H
