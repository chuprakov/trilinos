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

#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "NLPInterfacePack_NLPBarrier.hpp"
#include "MoochoPack_PostProcessBarrierLineSearch_Step.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )

namespace MoochoPack {

PostProcessBarrierLineSearch_Step::PostProcessBarrierLineSearch_Step(
  Teuchos::RefCountPtr<NLPInterfacePack::NLPBarrier> barrier_nlp
  )
	:
	barrier_nlp_(barrier_nlp)
	{
	TEST_FOR_EXCEPTION(
	  !barrier_nlp_.get(),
	  std::logic_error,
	  "PostProcessBarrierLineSearch_Step given NULL NLPBarrier."
	  );
	}
	

bool PostProcessBarrierLineSearch_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using Teuchos::dyn_cast;
	using IterationPack::print_algorithm_step;
	using AbstractLinAlgPack::Vp_StV;

	NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
	IpState             &s      = dyn_cast<IpState>(_algo.state());

   	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using IterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	// Get iteration quantities...
	value_type& f_kp1 = s.f().set_k(+1);
	f_kp1 = barrier_nlp_->objective_term();

	VectorMutable& Gf_kp1 = s.Gf().set_k(+1);
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
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	//const NLPAlgo   &algo = rsqp_algo(_algo);
	//const NLPAlgoState  &s    = algo.rsqp_state();
	out << L << "# Process out the temporary barrier term used for line search\n"
		<< L << "f_k -= barrier_term_k\n";
	}
} // end namespace MoochoPack 
