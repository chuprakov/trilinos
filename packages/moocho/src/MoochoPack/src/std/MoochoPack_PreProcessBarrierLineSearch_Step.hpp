// ////////////////////////////////////////////////////////////////////////////
// PreProcessBarrierLineSearch_Step.hpp
//
// Copyright (C) 2001
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

#ifndef PRE_PROCESS_BARRIER_LINE_SEARCH_STEP_H
#define PRE_PROCESS_BARRIER_LINE_SEARCH_STEP_H

#include <list>

#include "MoochoPack/src/MoochoPackTypes.hpp"
#include "MoochoPack/src/std/LineSearchFilter_Step.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "MoochoMoreUtilities/src/SetOptionsFromStreamNode.hpp"
#include "MoochoMoreUtilities/src/SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Fraction to boundary rule for calculating alpha max
 *   
 *
 * This class updates alpha_vl, alpha_vu, and alpha
 *  and x_kp1, Vl_kp1, and Vu_kp1
 *
 */

class PreProcessBarrierLineSearch_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		/** @name Constructors / initializers */
		//@{

		///
		/** Fraction to Boundary parameter
		 *
		 * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tau_boundary_frac )

		///
		/** Constructor.
		 */
		PreProcessBarrierLineSearch_Step(
		  Teuchos::RefCountPtr<NLPInterfacePack::NLPBarrier> barrier_nlp,
		  const value_type tau_boundary_frac = 0.99
		  );
		//@}

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
					 , poss_type assoc_step_poss);
		
		
		void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
						 , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
		//@}

	private:

    // //////////////////////////
    // Private data members

		Teuchos::RefCountPtr<NLPInterfacePack::NLPBarrier> barrier_nlp_;
		CastIQMember< Filter_T > filter_;

	}; // end class PreProcessBarrierLineSearch_Step


class PreProcessBarrierLineSearch_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode,
	  public OptionsFromStreamPack::SetOptionsToTargetBase< PreProcessBarrierLineSearch_Step >
	{
	public:
		PreProcessBarrierLineSearch_StepSetOptions(
		  PreProcessBarrierLineSearch_Step* target = 0,
		  const char opt_grp_name[] = "PreProcessBarrierLineSearch" );

	protected:

		/// Overridden from SetOptionsFromStreamNode
		void setOption( int option_num, const std::string& option_value );
	
	};	// end class PreProcessBarrierLineSearch_StepSetOptions

} // end namespace MoochoPack

#endif // #if !defined PRE_PROCESS_BARRIER_LINE_SEARCH_STEP_H
