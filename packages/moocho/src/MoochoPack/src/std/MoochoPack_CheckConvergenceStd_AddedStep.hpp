// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStep.hpp
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

#ifndef CHECK_CONVERGENCE_STD_ADDEDSTEP_H
#define CHECK_CONVERGENCE_STD_ADDEDSTEP_H

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.hpp"
#include "GeneralIterationPack/src/AlgorithmStep.hpp"
#include "StandardCompositionMacros.hpp"
#include "CheckConvergence_Strategy.hpp"

namespace ReducedSpaceSQPPack {

///
/** Check for convergence.
  */
class CheckConvergenceStd_AddedStep
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
	{
	public:
		
		///
		/** Strategy object to be used when checking for convergence
		 * 
		 */
		STANDARD_COMPOSITION_MEMBERS( CheckConvergence_Strategy, convergence_strategy )

		///
		CheckConvergenceStd_AddedStep(
		  MemMngPack::ref_count_ptr<CheckConvergence_Strategy> convergence_strategy
		  );

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
					 , poss_type assoc_step_poss);
		///
		void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
						 , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
		//@}

private:

};	// end class CheckConvergenceStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// CHECK_CONVERGENCE_STD_ADDEDSTEP_H
