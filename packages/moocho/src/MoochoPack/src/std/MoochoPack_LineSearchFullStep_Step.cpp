// ////////////////////////////////////////////////////////////////////////////
// LineSearchFullStep_Step.cpp
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

#include <ostream>

#include "ReducedSpaceSQPPack/src/std/LineSearchFullStep_Step.hpp"
#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackExceptions.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "GeneralIterationPack/src/print_algorithm_step.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpOut.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpOut.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"
#include "AbstractLinAlgPack/src/assert_print_nan_inf.hpp"
#include "ThrowException.hpp"

namespace ReducedSpaceSQPPack {

ReducedSpaceSQPPack::LineSearchFullStep_Step::LineSearchFullStep_Step(
		const bounds_tester_ptr_t&	bounds_tester
		)
	:
		bounds_tester_(bounds_tester)
{}


bool LineSearchFullStep_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_VpV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLP			&nlp	= algo.nlp();

	const size_type
		m = nlp.m();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}
	
	// alpha_k = 1.0
	IterQuantityAccess<value_type>
		&alpha_iq   = s.alpha();
	if( !alpha_iq.updated_k(0) )
		alpha_iq.set_k(0) = 1.0;

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nf_k        = " << s.f().get_k(0);
		if(m)
			out << "\n||c_k||inf = " << s.c().get_k(0).norm_inf();
		out << "\nalpha_k    = " << alpha_iq.get_k(0) << std::endl;
	}

	// x_kp1 = x_k + d_k
	IterQuantityAccess<VectorWithOpMutable>  &x_iq = s.x();
	const VectorWithOp                       &x_k   = x_iq.get_k(0);
	VectorWithOpMutable                      &x_kp1 = x_iq.set_k(+1);
	x_kp1 = x_k;
	Vp_StV( &x_kp1, alpha_iq.get_k(0), s.d().get_k(0) );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||x_kp1||inf   = " << s.x().get_k(+1).norm_inf() << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nx_kp1 =\n" << s.x().get_k(+1);
	}

	if(algo.algo_cntr().check_results()) {
		assert_print_nan_inf(
			x_kp1, "x_kp1",true
			,int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
			);
		if( nlp.num_bounded_x() ) {
			if(!bounds_tester().check_in_bounds(
				  int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
				, int(olevel) >= int(PRINT_VECTORS)					// print_all_warnings
				, int(olevel) >= int(PRINT_ITERATION_QUANTITIES)	// print_vectors
				, nlp.xl(), "xl"
				, nlp.xu(), "xu"
				, x_kp1, "x_kp1"
				))
			{
				THROW_EXCEPTION(
					true, TestFailed
					,"LineSearchFullStep_Step::do_step(...) : Error, "
					"the variables bounds xl <= x_k(+1) <= xu where violated!" );
			}
		}
	}

	// Calcuate f and c at the new point.
	nlp.set_multi_calc(true);
	nlp.set_f( &s.f().set_k(+1) );
	if(m) nlp.set_c( &s.c().set_k(+1) );
	nlp.calc_f(x_kp1);
	if(m) nlp.calc_c( x_kp1, false );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nf_kp1        = " << s.f().get_k(+1);
		if(m)
			out << "\n||c_kp1||inf = " << s.c().get_k(+1).norm_inf() << std::endl;
	}

	if( m && static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nc_kp1 =\n" << s.c().get_k(+1); 
	}

	if(algo.algo_cntr().check_results()) {
		assert_print_nan_inf( s.f().get_k(+1), "f(x_kp1)", true, &out );
		if(m)
			assert_print_nan_inf( s.c().get_k(+1), "c(x_kp1)", true, &out );
	}

	return true;
}

void LineSearchFullStep_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if alpha_k is not updated then\n"
		<< L << "    alpha_k = 1.0\n"
		<< L << "end\n"
		<< L << "x_kp1 = x_k + alpha_k * d_k\n"
		<< L << "f_kp1 = f(x_kp1)\n"
		<< L << "if m > 0 then c_kp1 = c(x_kp1)\n";
}

} // end namespace ReducedSpaceSQPPack
