// ////////////////////////////////////////////////////////////////////////////
// PreEvalNewPointBarrier_Step.cpp
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

#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/std/PreEvalNewPointBarrier_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"

#include "StringToBool.h"

#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

		/** Constructor.
		 */
PreEvalNewPointBarrier_Step::PreEvalNewPointBarrier_Step(
  const value_type relative_bound_push,
  const value_type absolute_bound_push
  )
	:
	relative_bound_push_(relative_bound_push),
	absolute_bound_push_(absolute_bound_push)
	{}

bool PreEvalNewPointBarrier_Step::do_step(
  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using DynamicCastHelperPack::dyn_cast;
	using GeneralIterationPack::print_algorithm_step;

	rSQPAlgo            &algo   = dyn_cast<rSQPAlgo>(_algo);
	ipState             &s      = dyn_cast<ipState>(_algo.state());
	NLPFirstOrderInfo   &nlp    = dyn_cast<NLPFirstOrderInfo>(algo.nlp());
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	if(!nlp.is_initialized())
		nlp.initialize(algo.algo_cntr().check_results());

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	IterQuantityAccess<VectorWithOpMutable> &x_iq = s.x();
	IterQuantityAccess<MatrixSymDiagonalStd> &Vl_iq = s.Vl();
	IterQuantityAccess<MatrixSymDiagonalStd> &Vu_iq = s.Vu();

	if( x_iq.last_updated() == IterQuantity::NONE_UPDATED ) 
		{
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
			{
			out << "\nInitialize x with x_k = nlp.xinit() ...\n"
				<< " and push x_k within bounds.\n";
			}
		VectorWithOpMutable& x_k = x_iq.set_k(0) = nlp.xinit();
  
		// apply transformation operator to push x sufficiently within bounds
		force_in_bounds_buffer(relative_bound_push_, 
							   absolute_bound_push_,
							   nlp.xl(),
							   nlp.xu(),
							   x_k);

		// evaluate the func and constraints
		IterQuantityAccess<VectorWithOpMutable>
			*c_iq   = nlp.m() > 0 ? &s.c() : NULL,
			*h_iq   = nlp.mI() > 0 ? &s.h() : NULL;

		using AbstractLinAlgPack::assert_print_nan_inf;
		if (assert_print_nan_inf(x_k, "x", true, NULL))
			{
			// Calcuate f and c at the new point.
			nlp.set_multi_calc(true);
			nlp.set_f( &s.f().set_k(0) );
			nlp.set_Gf( &s.Gf().set_k(0) );
			if (c_iq)
				{
				nlp.set_c( &c_iq->set_k(0) );
				nlp.set_Gc( &s.Gc().set_k(0) );
				}
			
			if (h_iq)
				{
				nlp.set_h( &h_iq->set_k(0) );
				nlp.set_Gh( &s.Gh().set_k(0) );
				}

			nlp.calc_Gf(x_k, true);
			nlp.calc_f( x_k, false ); 

			if (c_iq)
				{
				nlp.calc_Gc( x_k, false );
				nlp.calc_c( x_k, false);
				}
			
			if (h_iq)
				{
				nlp.calc_Gh( x_k, false );
				nlp.calc_h( x_k, false );
				}
			}
		}

	if (s.barrier_parameter().last_updated() == IterQuantity::NONE_UPDATED)
		{
		s.barrier_parameter().set_k(-1) = 0.1;
		}

	// Print vector information
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "x_k=\n";
		x_iq.get_k(0).output(out);
		}

	return true;
	}


void PreEvalNewPointBarrier_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	//const rSQPAlgo   &algo = rsqp_algo(_algo);
	//const rSQPState  &s    = algo.rsqp_state();
	out << L << "# Evaluate information specific to primal / dual barrier algorithms\n"
		<< L << "if (x never updated) then\n"
		<< L << "   x_k = nlp.xinit()\n"
		<< L << "   force_in_bounds(x_k)\n"
		<< L << "   Gf_k = calc_Gf\n"
		<< L << "   f_k = calc_f\n"
		<< L << "   Gc_k = calc_Gf\n"
		<< L << "   c_k = calc_f\n"
		<< L << "   Gh_k = calc_Gf\n"
		<< L << "   h_k = calc_f\n"
		<< L << "end\n";
	}


namespace {

const int local_num_options = 2;

enum local_EOptions 
	{
	RELATIVE_BOUND_PUSH,
	ABSOLUTE_BOUND_PUSH
	};

const char* local_SOptions[local_num_options] = 
	{
	"relative_bound_push",
	"absolute_bound_push"
	};

}

 
PreEvalNewPointBarrier_StepSetOptions::PreEvalNewPointBarrier_StepSetOptions(
  PreEvalNewPointBarrier_Step* target
  , const char opt_grp_name[] )
	:
	OptionsFromStreamPack::SetOptionsFromStreamNode(
	  opt_grp_name, local_num_options, local_SOptions ),
	OptionsFromStreamPack::SetOptionsToTargetBase< PreEvalNewPointBarrier_Step >( target )
	{
	}

void PreEvalNewPointBarrier_StepSetOptions::set_option( 
  int option_num, const std::string& option_value )
	{
	using OptionsFromStreamPack::StringToBool;

	typedef PreEvalNewPointBarrier_Step target_t;
	switch( (local_EOptions)option_num ) 
		{
		case RELATIVE_BOUND_PUSH:
			target().relative_bound_push(::atof(option_value.c_str()));
			break;
		case ABSOLUTE_BOUND_PUSH:
			target().absolute_bound_push(::atof(option_value.c_str()));
			break;
		default:
			assert(0);	// Local error only?
		}
	}

}; // end namespace ReducedSpaceSQPPack
