// ////////////////////////////////////////////////////////////////////////////
// PostProcessBarrierLineSearch_Step.cpp
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

#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "NLPInterfacePack/include/BarrierNLP.h"
#include "ReducedSpaceSQPPack/include/std/PostProcessBarrierLineSearch_Step.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )

namespace ReducedSpaceSQPPack {

PostProcessBarrierLineSearch_Step::PostProcessBarrierLineSearch_Step(
  MemMngPack::ref_count_ptr<NLPInterfacePack::BarrierNLP> barrier_nlp
  )
	:
	barrier_nlp_(barrier_nlp)
	{
	THROW_EXCEPTION(
	  !barrier_nlp_.get(),
	  std::logic_error,
	  "PostProcessBarrierLineSearch_Step given NULL BarrierNLP."
	  );
	}
	

bool PostProcessBarrierLineSearch_Step::do_step(
  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using DynamicCastHelperPack::dyn_cast;
	using GeneralIterationPack::print_algorithm_step;
	using AbstractLinAlgPack::Vp_StV;

	rSQPAlgo            &algo   = dyn_cast<rSQPAlgo>(_algo);
	ipState             &s      = dyn_cast<ipState>(_algo.state());

   	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	// Get iteration quantities...
	value_type& f_kp1 = s.f().set_k(+1);
	f_kp1 = barrier_nlp_->objective_term();

	VectorWithOpMutable& Gf_kp1 = s.Gf().set_k(+1);
	Gf_kp1 = *(barrier_nlp_->grad_objective_term());
	
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "\nf = " << f_kp1
			<< "\nbarrier_term = " << barrier_nlp_->barrier_term() << std::endl;
		
		}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "Gf = \n" << Gf_kp1
			<< "\ngrad_barrier_term  = \n" << *(barrier_nlp_->grad_barrier_term());
		
		}
	return true;
	}


void PostProcessBarrierLineSearch_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	//const rSQPAlgo   &algo = rsqp_algo(_algo);
	//const rSQPState  &s    = algo.rsqp_state();
	out << L << "# Process out the temporary barrier term used for line search\n"
		<< L << "f_k -= barrier_term_k\n";
	}
} // end namespace ReducedSpaceSQPPack 
