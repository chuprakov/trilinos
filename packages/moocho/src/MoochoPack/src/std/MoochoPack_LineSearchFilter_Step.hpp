// ////////////////////////////////////////////////////////////////////////////
// LineSearchFilter_Step.hpp
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

#ifndef LINE_SEARCH_FILTER_STEP_H
#define LINE_SEARCH_FILTER_STEP_H

#include <list>

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.hpp"
#include "GeneralIterationPack/src/AlgorithmStep.hpp"
#include "GeneralIterationPack/src/CastIQMember.hpp"
#include "StandardCompositionMacros.hpp"
#include "StandardMemberCompositionMacros.hpp"

#include "GeneralIterationPack/src/AlgorithmState.hpp"

namespace ReducedSpaceSQPPack {
 
// structure for storing filter points 
class FilterEntry 
	{
	public:
		FilterEntry( value_type new_f, value_type new_theta, int new_iter)
			: f(new_f), theta(new_theta), iter(new_iter) {}
		
		value_type f;
		value_type theta;
		int iter;
	};

typedef std::list< FilterEntry > Filter_T;
const std::string FILTER_IQ_STRING = "LS_FilterEntries";


///
/** Filter line-search step class.
 *
 * Todo: Finish documentataion.
 */
class LineSearchFilter_Step
    : public GeneralIterationPack::AlgorithmStep // doxygen needs full path
	{
	public:
		
		/** @name Public types */
		//@{
	
		//@}
	
		/** @name Constructors / initializers */
		//@{
	
		///
		/** Feasibility decrease fraction.
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, gamma_theta )

		///
		/** Optimality decrease fraction.
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, gamma_f )

		///
		/** alpha_min linearization correction fraction
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, gamma_alpha )

		///
		/** Delta parameter for switching condition.
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, delta )
	
		///
		/** Exponent for objective in switching condition.
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, s_f )

		///
		/** Exponent for theta in switching condition.
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, s_theta )

		///
		/** Factor to evaluate theta_small
		 * theta_small = theta_small_fact*max(1,theta_k)
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_small_fact )

		///
		/** Constant for Armijo condition on objective
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta_f )

		///
		/** Backtracking fraction for step.
		 *
		 * ToDo: Finish documentation.
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, back_track_frac )

		///
		/** Constructor.
		 */
		LineSearchFilter_Step(
		  MemMngPack::ref_count_ptr<NLPInterfacePack::NLP> nlp
		  ,const std::string obj_iq_name = "f"
		  ,const std::string grad_obj_iq_name = "Gf"
		  ,const value_type           &gamma_theta      = 1e-5
		  ,const value_type          &gamma_f          = 1e-5
		  ,const value_type          &gamma_alpha      = 5e-2
		  ,const value_type          &delta            = 1e-4
		  ,const value_type          &s_theta          = 1.1
		  ,const value_type          &s_f              = 2.3
		  ,const value_type          &theta_small_fact = 1e-4 
		  ,const value_type          &eta_f            = 1e-4
		  ,const value_type          &back_track_frac  = 0.5
		  );

		//@}
	
		///
		/** Destructor.
		 */
		~LineSearchFilter_Step();

		//@}

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step( Algorithm& algo, poss_type step_poss,
					  GeneralIterationPack::EDoStepType type,
					  poss_type assoc_step_poss);
		///
		void print_step( const Algorithm& algo, poss_type step_poss,
						 GeneralIterationPack::EDoStepType type,
						 poss_type assoc_step_poss, std::ostream& out,
						 const std::string& leading_str ) const;
		//@}

	private:

		// Not defined and not to be called
		//	LineSearchFilter_Step();

		// Private Data
		CastIQMember<Filter_T> filter_;

		/// Iteration quantity access for objective value
		CastIQMember<value_type> obj_f_;

		/// ITeration quantity access for objective gradient
		CastIQMember<VectorWithOpMutable> grad_obj_f_;

		

		// nlp to use for calculations
		MemMngPack::ref_count_ptr<NLPInterfacePack::NLP> nlp_;

		// Validate input parameters - fix if possible
		bool ValidatePoint( IterQuantityAccess<VectorWithOpMutable>& x,
							IterQuantityAccess<value_type>& f,
							IterQuantityAccess<VectorWithOpMutable>* c,
							IterQuantityAccess<VectorWithOpMutable>* h,
							bool throw_excpt ) const;


		// Check that new point is not within taboo region of filter
		bool CheckFilterAcceptability( value_type f, 
									   value_type theta,
									   AlgorithmState& s ) const;


		// Check the Armijo condition on f
		bool CheckArmijo( value_type Gf_t_dk, 
						  value_type alpha_k, 
						  const IterQuantityAccess<value_type>& f_iq ) const;


		// Check if f or c has sufficient reduction
		bool CheckFractionalReduction( const IterQuantityAccess<value_type>& f_iq,
									   value_type theta_kp1, 
									   value_type theta_k ) const;


		// Calculate the new point given d and alpha
		void UpdatePoint( const VectorWithOpMutable& d,
						  value_type alpha, 
						  IterQuantityAccess<VectorWithOpMutable>& x,
						  IterQuantityAccess<value_type>& f,
						  IterQuantityAccess<VectorWithOpMutable>* c,
						  IterQuantityAccess<VectorWithOpMutable>* h,
						  NLP& nlp ) const;


		// Calculate the minimum alpha before hitting restoration phase
		value_type CalculateAlphaMin( value_type Gf_t_dk, 
									  value_type theta_k,
									  value_type theta_small ) const;


		// Calculate the constraint norm
		value_type CalculateTheta_k( IterQuantityAccess<VectorWithOpMutable>* c,
									 IterQuantityAccess<VectorWithOpMutable>* h,
									 int k ) const;


		// decide if we should switch to Armijo for objective
		bool ShouldSwitchToArmijo( value_type Gf_t_dk,
								   value_type alpha_k,
								   value_type theta_k,
								   value_type theta_small) const;

		// Update the filter from the last iteration
		void UpdateFilter( GeneralIterationPack::AlgorithmState& s ) const;

		// Update the filter from the last iteration and Augment it with
		// the new point
		void AugmentFilter( value_type f_with_boundary,
							value_type theta_with_boundary,
							GeneralIterationPack::AlgorithmState& s ) const;

	};	// end class LineSearchFilter_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// LINE_SEARCH_FILTER_STEP_H
