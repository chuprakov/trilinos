// ////////////////////////////////////////////////////////////////////////////
// PostEvalNewPointBarrier_Step.h
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

#ifndef POST_EVAL_NEW_POINT_BARRIER_STEP_H
#define POST_EVAL_NEW_POINT_BARRIER_STEP_H


#include "GeneralIterationPack/src/AlgorithmStep.h"
#include "StandardCompositionMacros.h"
#include "StandardMemberCompositionMacros.h"

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/src/AlgorithmStep.h"

#include "SetOptionsFromStreamNode.h"
#include "SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Standard evaluation step class for extra parameters in primal/dual barrier method.
 *
 * This class calculates \c invXu, \c invXl \c invXu_m_invXl
 *
 */

class PostEvalNewPointBarrier_Step
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
					 , poss_type assoc_step_poss);
		
		
		void print_step( const GeneralIterationPack::Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
						 , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
		//@}


	}; // end class PostEvalNewPointBarrier_Step

}  // end namespace ReducedSpaceSQPPack

#endif // POST_EVAL_NEW_POINT_BARRIER_STEP_H
