//////////////////////////////////////////////////////////////////////
// MoochoPack_UpdateReducedSigma_Step.hpp
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

#ifndef UPDATE_REDUCED_SIGMA_STEP_H
#define UPDATE_REDUCED_SIGMA_STEP_H

#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Standard class for updating the reduced sigma 
 *   for interior point optimization
 *
 *
 */

class UpdateReducedSigma_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		enum e_update_methods
			{
			ALWAYS_EXPLICIT,
			BFGS_PRIMAL,
			BFGS_DUAL_NO_CORRECTION,
			BFGS_DUAL_EXPLICIT_CORRECTION,
			BFGS_DUAL_SCALING_CORRECTION
			};

		///
		/** update method for the reduced sigma term
		 *  update_method = always_explicit;
		 *	update_method = BFGS_primal;
		 *	update_method = BFGS_dual_no_correction;
		 *	update_method = BFGS_dual_explicit_correction; *** (default)
		 *	update_method = BFGS_dual_scaling_correction;
		 *** These options determine exactly how the reduced sigma
		 *** term will be updated. 
		 ***
		 *** always_explicit 				: the full Z_kT*Sigma*Zk at each step (expensive)
		 *** BFGS_primal 					: a BFGS update of mu*X^-2 (exact at solution)
		 *** BFGS_dual_no_correction 		: update with Z_kT*Sigma*Z_k*pz 
		 ***										(no correction when mu changes)
		 *** BFGS_dual_explicit_correction 	: same as above
		 ***										(do an explicit calculation when mu changes)
		 *** BFGS_dual_scaling_correction    : same as above
		 ***										(scale by mu_kp1/mu_k when mu changes)
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( e_update_methods, update_method)

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
		UpdateReducedSigma_Step(
//		  const e_update_methods update_method = BFGS_DUAL_EXPLICIT_CORRECTION
		  const e_update_methods update_method = ALWAYS_EXPLICIT // For now only!
		  );
		//@}

	private:
		void FormReducedSigmaExplicitly(NLPAlgo& algo, IpState& s, EJournalOutputLevel olevel,  std::ostream& out);


	}; // end class EvalNewPointBarrier_Step

const char UpdateReducedSigma_opt_grp_name[] = "UpdateReducedSigma";
class UpdateReducedSigma_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode,
	  public OptionsFromStreamPack::SetOptionsToTargetBase< UpdateReducedSigma_Step >
	{
	public:
		UpdateReducedSigma_StepSetOptions(
		  UpdateReducedSigma_Step* target = 0,
		  const char opt_grp_name[] = UpdateReducedSigma_opt_grp_name );

	protected:

		/// Overridden from SetOptionsFromStreamNode
		void setOption( int option_num, const std::string& option_value );
	
	};	// end class UpdateReducedSigma_StepSetOptions


}  // end namespace MoochoPack

#endif // UPDATE_REDUCED_SIGMA_STEP_H
