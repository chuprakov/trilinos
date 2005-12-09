// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_UpdateBarrierParameter_Step.hpp
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

#ifndef UPDATE_BARRIER_PARAMETER_STEP_H
#define UPDATE_BARRIER_PARAMETER_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Barrier Parameter (mu) Update 
 *
 * This class updates barrier_parameter & e_tol for next iteration
 *
 */

class UpdateBarrierParameter_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		/** @name Constructors / initializers */
		//@{

		///
		/** Initial barrier parameter
		 *
		 * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, init_barrier_parameter )

		///
		/** barrier_parameter decrease fraction (linear decrease)
		 *
		 * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tau_mu )

		///
		/** barrier_parameter decrease power (for superlinear decrease)
		 *
		 * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_mu )

		///
		/** error tolerance fraction
		 *
		 * e_tol = min( e_tol_max, tau_epsilon*min(mu_k,mu_k^theta_epsilon))
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tau_epsilon )

		///
		/** error tolerance power
		 *
		 * e_tol = min( e_tol_max, tau_epsilon*min(mu_k,mu_k^theta_epsilon))
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_epsilon )

		///
		/** maximum error tolerance
		 *
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, e_tol_max )

		///
		/** Constructor.
		 */
		UpdateBarrierParameter_Step(
		  const value_type init_barrier_parameter = 0.1,
		  const value_type tau_mu = 0.2,
		  const value_type theta_mu = 1.5,
		  const value_type tau_epsilon = 10,
		  const value_type theta_epsilon = 1.1,
		  const value_type e_tol_max = 1000
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
		value_type Calculate_e_tol(value_type mu);

	}; // end class UpdateBarrierParameter_Step


class UpdateBarrierParameter_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode,
	  public OptionsFromStreamPack::SetOptionsToTargetBase< UpdateBarrierParameter_Step >
	{
	public:
		UpdateBarrierParameter_StepSetOptions(
		  UpdateBarrierParameter_Step* target = 0,
		  const char opt_grp_name[] = "UpdateBarrierParameter" );

	protected:

		/// Overridden from SetOptionsFromStreamNode
		void setOption( int option_num, const std::string& option_value );
	
	};	// end class UpdateBarrierParameter_StepSetOptions

} // end namespace MoochoPack

#endif // #if !defined UPDATE_BARRIER_PARAMETER_STEP_H
