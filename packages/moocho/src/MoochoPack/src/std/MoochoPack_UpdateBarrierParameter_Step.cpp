// ////////////////////////////////////////////////////////////////////////////
// UpdateBarrierParameter_Step.cpp
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

#include <ostream>
#include <typeinfo>
#include <iostream>
#include <math.h>

#include "ReducedSpaceSQPPack/src/std/UpdateBarrierParameter_Step.h"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.h"
#include "GeneralIterationPack/src/print_algorithm_step.h"
#include "dynamic_cast_verbose.h"
#include "ReducedSpaceSQPPack/src/ipState.h"
#include "AbstractLinAlgPack/src/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/src/VectorStdOps.h"
#include "AbstractLinAlgPack/src/VectorAuxiliaryOps.h"
#include "ThrowException.h"

#include "StringToBool.h"

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )

namespace ReducedSpaceSQPPack {

UpdateBarrierParameter_Step::UpdateBarrierParameter_Step(
  const value_type init_barrier_parameter,
  const value_type tau_mu,
  const value_type theta_mu,
  const value_type tau_epsilon,
  const value_type theta_epsilon,
  const value_type e_tol_max
  )
	:
	init_barrier_parameter_(init_barrier_parameter),
	tau_mu_(tau_mu),
	theta_mu_(theta_mu),
	tau_epsilon_(tau_epsilon),
	theta_epsilon_(theta_epsilon),
	e_tol_max_(e_tol_max)
	{}
	

bool UpdateBarrierParameter_Step::do_step(
  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using DynamicCastHelperPack::dyn_cast;
	using GeneralIterationPack::print_algorithm_step;

	rSQPAlgo            &algo   = dyn_cast<rSQPAlgo>(_algo);
	ipState             &s      = dyn_cast<ipState>(_algo.state());
	NLP                 &nlp    = algo.nlp();

	if (!nlp.is_initialized())
		{
		nlp.initialize(false);
		}
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	
	///***********************************************************
	// Get iteration quantities
	///***********************************************************
	IterQuantityAccess<value_type>& e_tol_iq = s.e_tol();
	IterQuantityAccess<value_type>& mu_iq = s.barrier_parameter();

	///***********************************************************
	// Check values and initialize, if necessary
	///***********************************************************
	/*	if (mu_iq.last_updated() == IterQuantity::NONE_UPDATED)
		{
		mu_iq.set_k(0) = init_barrier_parameter_;
		e_tol_iq.set_k(0) = Calculate_e_tol(mu_iq.get_k(0));

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
			{
			out << "\nInitializing barrier parameter (mu) and sub problem tolerance (e_tol) ...\n";
			}
		}
		else*/
		{
		///***********************************************************
		// if tolerance is satisfied, calculate new barrier_parameter 
		//  and e_tol, otherwise update mu and e_tol from last iter
		///***********************************************************
		const value_type opt_err = s.opt_kkt_err().get_k(0);
		const value_type feas_err = s.feas_kkt_err().get_k(0);
		const value_type comp_err = s.comp_kkt_err().get_k(0);

		const value_type mu_km1 = mu_iq.get_k(-1);
		if (e_tol_iq.last_updated() == IterQuantity::NONE_UPDATED)
			{
			// First time through, don't let mu change
			mu_iq.set_k(0,-1);
			e_tol_iq.set_k(0) = Calculate_e_tol(mu_iq.get_k(0));
        	}
		else
			{
			const value_type e_tol_km1 = e_tol_iq.get_k(-1);
			bool sub_prob_converged = (opt_err < e_tol_km1 && feas_err < e_tol_km1 && comp_err < e_tol_km1); 
			if (sub_prob_converged)
				{
				// Calculate new mu and e_tol
				value_type& mu_k = mu_iq.set_k(0);
				mu_k = min(tau_mu_*mu_km1, pow(mu_km1, theta_mu_));
				e_tol_iq.set_k(0) = Calculate_e_tol(mu_k);
				
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
					{
					out << "\nSub-problem converged!\n"
						<< " Updating barrier parameter (mu) and sub problem tolerance (e_tol) ...\n";
					}


	        	        /*VectorWithOpMutable& vu = s.Vu().set_k(0).diag();
        	        	vu = 0;
		                Vp_StV(&vu, mu_k, s.invXu().get_k(0).diag());
                		correct_upper_bound_multipliers(nlp.xu(), NLP::infinite_bound(), &vu); 
                
		                VectorWithOpMutable& vl = s.Vl().set_k(0).diag();
                		vl = 0;
		                Vp_StV(&vl, mu_k, s.invXl().get_k(0).diag());
                		correct_lower_bound_multipliers(nlp.xl(), -NLP::infinite_bound(), &vl);*/
				}
			else
				{
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
					{
					out << "\nSub-problem not-converged!\n"
						<< " Keeping existing barrier parameter (mu) and sub problem tolerance (e_tol) ...\n";
					}
				mu_iq.set_k(0,-1);
				e_tol_iq.set_k(0,-1);
				}
			}

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
			{
			out << "\nbarrier_parameter (mu) = " << mu_iq.get_k(0)
				<< "\nsub problem tolerance (e_tol) = " << e_tol_iq.get_k(0)  << std::endl;
			}
		}

	return true;
	}


void UpdateBarrierParameter_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	//const rSQPAlgo   &algo = rsqp_algo(_algo);
	//const rSQPState  &s    = algo.rsqp_state();
	out << L << "# Update the interior point barrier parameter (mu)\n"
		<< L << "if (KKTerror < e_tol) then\n"
		<< L << "   mu_kp1 = min(tau_mu*mu_k, mu_k^theta_mu)\n"
		<< L << "   e_tol_kp1 = min(e_tol_max, tau_epsilon*min(mu_k, mu_k^theta_epsilon))\n"
		<< L << "else\n"
		<< L << "   mu_kp1 = mu_k\n"
		<< L << "   e_tol_kp1 = e_tol_k\n"
		<< L << "end;\n";
	}

value_type UpdateBarrierParameter_Step::Calculate_e_tol(value_type mu)
	{	
	value_type e_tol = tau_epsilon_*min(mu, pow(mu, theta_epsilon_));
	e_tol = min(e_tol_max_, e_tol);

	return e_tol;
	}


namespace {

const int local_num_options = 5;

enum local_EOptions 
	{
		TAU_MU,
		THETA_MU,
		TAU_EPSILON,
		THETA_EPSILON,
		E_TOL_MAX
	};

const char* local_SOptions[local_num_options] = 
	{
		"tau_mu",
		"theta_mu",
		"tau_epsilon",
		"theta_epsilon",
		"e_tol_max"
	};

}

 
UpdateBarrierParameter_StepSetOptions::UpdateBarrierParameter_StepSetOptions(
  UpdateBarrierParameter_Step* target
  , const char opt_grp_name[] )
	:
	OptionsFromStreamPack::SetOptionsFromStreamNode(
	  opt_grp_name, local_num_options, local_SOptions ),
	OptionsFromStreamPack::SetOptionsToTargetBase< UpdateBarrierParameter_Step >( target )
	{
	}

void UpdateBarrierParameter_StepSetOptions::set_option( 
  int option_num, const std::string& option_value )
	{
	using OptionsFromStreamPack::StringToBool;
  
	typedef UpdateBarrierParameter_Step target_t;
	switch( (local_EOptions)option_num ) 
		{
		case TAU_MU:
			target().tau_mu(::atof(option_value.c_str()));
			break;
		case THETA_MU:
			target().theta_mu(::atof(option_value.c_str()));
			break;
		case TAU_EPSILON:
			target().tau_epsilon(::atof(option_value.c_str()));
			break;
		case THETA_EPSILON:
			target().theta_epsilon(::atof(option_value.c_str()));
			break;
		case E_TOL_MAX:
			target().e_tol_max(::atof(option_value.c_str()));
			break;
		default:
			assert(0);	// Local error only?
		}
	}

} // end namespace ReducedSpaceSQPPack 
