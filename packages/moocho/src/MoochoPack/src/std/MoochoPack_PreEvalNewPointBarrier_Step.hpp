// ////////////////////////////////////////////////////////////////////////////
// PreEvalNewPointBarrier_Step.hpp
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

#ifndef PRE_EVAL_NEW_POINT_BARRIER_STEP_H
#define PRE_EVAL_NEW_POINT_BARRIER_STEP_H


#include "IterationPack/src/AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "MoochoPack/src/MoochoPackTypes.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"

#include "MoochoMoreUtilities/src/SetOptionsFromStreamNode.hpp"
#include "MoochoMoreUtilities/src/SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Standard evaluation step class for extra parameters in primal/dual barrier method.
 *
 * This class calculates \c invXu, \c invXl \c invXu_m_invXl
 *
 */

class PreEvalNewPointBarrier_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		///
		/** relative fraction for initializing x within
		 *   bounds.
		 *   xl_sb = min(xl+relative_bound_push*(xu-xl),
		 *               xl + absolute_bound_push)
		 *   xu_sb = max(xu-relative_bound_push*(xu-xl),
		 *               xu - absolute_bound_push)
		 *   if (xl_sb > xu_sb) then
		 *      x = (xl + (xu-xl)/2
		 *   else if (x < xl_sb) then 
		 *      x = xl_sb
		 *   else if (x > xu_sb) then
		 *      x = xu_sb
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, relative_bound_push )

		///
		/** absolute measure for initializing x within
		 *   bounds.
		 *   xl_sb = min(xl+relative_bound_push*(xu-xl),
		 *               xl + absolute_bound_push)
		 *   xu_sb = max(xu-relative_bound_push*(xu-xl),
		 *               xu - absolute_bound_push)
		 *   if (xl_sb > xu_sb) then
		 *      x = (xl + (xu-xl)/2
		 *   else if (x < xl_sb) then 
		 *      x = xl_sb
		 *   else if (x > xu_sb) then
		 *      x = xu_sb
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, absolute_bound_push )

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
					 , poss_type assoc_step_poss);
		
		
		void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
						 , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
		//@}

		/** Constructor.
		 */
		PreEvalNewPointBarrier_Step(
		  const value_type relative_bound_push = 0.01,
		  const value_type absolute_bound_push = 0.001
		  );
		//@}

	}; // end class PreEvalNewPointBarrier_Step

class PreEvalNewPointBarrier_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode,
	  public OptionsFromStreamPack::SetOptionsToTargetBase< PreEvalNewPointBarrier_Step >
	{
	public:
		PreEvalNewPointBarrier_StepSetOptions(
		  PreEvalNewPointBarrier_Step* target = 0,
		  const char opt_grp_name[] = "PreEvalNewPointBarrier" );

	protected:

		/// Overridden from SetOptionsFromStreamNode
		void setOption( int option_num, const std::string& option_value );
	
	};	// end class PreEvalNewPointBarrier_StepSetOptions


}  // end namespace MoochoPack

#endif // PRE_EVAL_NEW_POINT_BARRIER_STEP_H
