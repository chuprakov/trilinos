// ////////////////////////////////////////////////////////////////////////////
// PostProcessBarrierLineSearch_Step.hpp
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

#ifndef POST_PROCESS_BARRIER_LINE_SEARCH_STEP_H
#define POST_PROCESS_BARRIER_LINE_SEARCH_STEP_H

#include "MoochoPack/src/MoochoPackTypes.hpp"
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

class PostProcessBarrierLineSearch_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		/** @name Constructors / initializers */
		//@{

		///
		/** Constructor.
		 */
		PostProcessBarrierLineSearch_Step(
		  Teuchos::RefCountPtr<NLPInterfacePack::NLPBarrier> barrier_nlp
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
		Teuchos::RefCountPtr<NLPInterfacePack::NLPBarrier> barrier_nlp_;

	}; // end class PostProcessBarrierLineSearch_Step


} // end namespace MoochoPack

#endif // #if !defined POST_PROCESS_BARRIER_LINE_SEARCH_STEP_H
