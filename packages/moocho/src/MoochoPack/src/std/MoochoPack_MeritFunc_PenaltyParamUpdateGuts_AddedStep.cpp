// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateGuts_AddedStep.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateGuts_AddedStep.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "ConstrainedOptimizationPack/include/MeritFuncPenaltyParam.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPDirecDeriv.h"
#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"

namespace ReducedSpaceSQPPack {

MeritFunc_PenaltyParamUpdateGuts_AddedStep::MeritFunc_PenaltyParamUpdateGuts_AddedStep(
	value_type                     small_mu
	,value_type                    mult_factor
	,value_type                    kkt_near_sol
	)
	:near_solution_(false)
	,small_mu_(small_mu)
	,mult_factor_(mult_factor)
	,kkt_near_sol_(kkt_near_sol)
{}

bool MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLP			&nlp	= algo.nlp();
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	const size_type
		n  = nlp.n(),
		m  = nlp.m(),
		mI = nlp.mI();
	IterQuantityAccess<MeritFuncNLP>
		&merit_func_nlp_iq = s.merit_func_nlp();

	if( !merit_func_nlp_iq.updated_k(0) ) {
		const int merit_func_k_last_updated = merit_func_nlp_iq.last_updated();
		if( merit_func_k_last_updated != IterQuantity::NONE_UPDATED ) {
			MeritFuncNLP
				&merit_func_nlp_k_last = merit_func_nlp_iq.get_k(merit_func_k_last_updated);
			merit_func_nlp_iq.set_k(0) = merit_func_nlp_k_last;
		}
		else {
			merit_func_nlp_iq.set_k(0); // Just use default constructor
		}
		MeritFuncNLP
			&merit_func_nlp_k = merit_func_nlp_iq.get_k(0);
		MeritFuncPenaltyParam
			*param = dynamic_cast<MeritFuncPenaltyParam*>(&merit_func_nlp_k);
		THROW_EXCEPTION(
			!param, std::logic_error
			,"MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(...), Error "
			<< "The class " << typeid(merit_func_nlp_k).name() << " does not support the "
			<< "MeritFuncPenaltyParam iterface" );
		MeritFuncNLPDirecDeriv
			*direc_deriv = dynamic_cast<MeritFuncNLPDirecDeriv*>(&merit_func_nlp_k);
		THROW_EXCEPTION(
			!direc_deriv, std::logic_error
			,"MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(...), Error "
			<< "The class " << typeid(merit_func_nlp_k).name() << " does not support the "
			<< "MeritFuncNLPDirecDeriv iterface" );
		value_type  new_mu = 0.0;
		value_type  min_mu = 0.0;
		if ( this->min_mu(s,&min_mu) ) {
			// Update the penalty parameter as defined in the fortran rSQP code (EXACT2())
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nUpdate the penalty parameter...\n";
			}
			value_type
				mu_km1 = param->mu(),
				mult_fact = (1.0 + mult_factor_);
			if(near_solution_) {
				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
					out << "\nNear solution, forcing mu_k >= mu_km1...\n";
				}
				new_mu = std::_MAX( std::_MAX( mu_km1, mult_fact * min_mu ), small_mu_ );
			}
			else {
				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
					out << "\nNot near solution, allowing reduction in mu ...\n";
				}
				new_mu =	std::_MAX(
					(3.0 * mu_km1 + min_mu) / 4.0	
					, std::_MAX( mult_fact * min_mu, small_mu_ )
					); 
				value_type
					kkt_error = s.opt_kkt_err().get_k(0) + s.feas_kkt_err().get_k(0);
				if(kkt_error <= kkt_near_sol_) {
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out << "\nkkt_error = " << kkt_error << " <= kkt_near_sol = "
							<< kkt_near_sol_ << std::endl
							<< "Switching to forcing mu_k >= mu_km1 in the future\n";
					}
					near_solution_ = true;
				}
			}
		}
		else {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nDon't have the info to update penalty parameter so just use the last updated...\n";
			}
			new_mu = param->mu();
		}
		// Set the penalty parameter
		param->mu( new_mu );
		// In addition also compute the directional derivative
		direc_deriv->calc_deriv(
			s.Gf().get_k(0)
			,m  ? &s.c().get_k(0) : NULL
			,mI ? &s.h().get_k(0) : NULL
			,mI ? &nlp.hl()       : NULL
			,mI ? &nlp.hu()       : NULL
			,s.d().get_k(0)
			);
		
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\nmu = " << new_mu << "\n";
		}
	}
	return true;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Update the penalty parameter for the merit function to ensure\n"
		<< L << "*** a descent direction a directional derivatieve.\n"
		<< L << "*** phi is a merit function object that uses the penalty parameter mu.\n"
		<< L << "default: near_solution = false\n"
		<< L << "         small_mu = " << small_mu_ << std::endl
		<< L << "         mult_factor = " << mult_factor_ << std::endl
		<< L << "         kkt_near_sol = " << kkt_near_sol_ << std::endl
		<< L << "if merit_func_nlp_k is not already updated then\n"
		<< L << "    if some merit_func_nlp_k(?) has been udpated then\n"
		<< L << "        merit_func_nlp_k = merit_func_nlp_k(last_udpated)\n"
		<< L << "    else\n"
		<< L << "        merit_func_nlp_k = default construction\n"
		<< L << "    end\n"
		<< L << "    if merit_func_nlp_k does not support MeritFuncPenaltyParam throw excpetion\n"
		<< L << "    if merit_func_nlp_k does not support MeritFuncNLPDirecDeriv throw excpetion\n"
		;
	            print_min_mu_step( out, L + "    " ); 
	out
		<< L << "   mu_new = merit_func_nlp_k.mu()\n"
		<< L << "   if update_mu == true then\n"
		<< L << "       mu_last = merit_func_nlp_k.mu()\n"
		<< L << "       mult_fact = 1.0 + mult_factor\n"
		<< L << "       if near_solution == true\n"
		<< L << "           mu_new = max( max( mu_last, mult_fact*min_mu ), small_mu )\n"
		<< L << "       else\n"
		<< L << "           mu_new = max(   ( 3.0 * mu_last + min_mu ) / 4.0\n"
		<< L << "                           , max( mult_fact * min_mu , small_mu ) )\n"
		<< L << "           kkt_error = opt_kkt_err_k + feas_kkt_err_k\n"
		<< L << "           if kkt_error <= kkt_near_sol then\n"
		<< L << "               near_solution = true\n"
		<< L << "           end\n"
		<< L << "       end\n"
		<< L << "   else\n"
		<< L << "       mu_new = merit_func_nlp_k.mu()\n"
		<< L << "   end\n"
		<< L << "   merit_func_nlp_k..mu(mu_new)\n"
		<< L << "   merit_func_nlp_k.calc_deriv(Gf_k,c_k,h_k,hl,hu,d_k)\n"
		<< L << "end\n"
		;
}

// Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::small_mu( value_type small_mu )
{
	small_mu_ = small_mu;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::small_mu() const
{
	return small_mu_;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::mult_factor( value_type mult_factor )
{
	mult_factor_ = mult_factor;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::mult_factor() const
{
	return mult_factor_;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::kkt_near_sol( value_type kkt_near_sol )
{
	kkt_near_sol_ = kkt_near_sol;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::kkt_near_sol() const
{
	return kkt_near_sol_;
}

}	// end namespace ReducedSpaceSQPPack
