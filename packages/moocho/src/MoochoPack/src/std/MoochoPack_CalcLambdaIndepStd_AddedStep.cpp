// ////////////////////////////////////////////////////////////////////////////
// CalcLambdaIndepStd_AddedStep.cpp
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

#include "../std/CalcLambdaIndepStd_AddedStep.h"
#include "../rsqp_algo_conversion.h"
#include "GeneralIterationPack/src/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/src/ComputeMinMult.h"
#include "ConstrainedOptimizationPack/src/VectorWithNorms.h"
#include "SparseLinAlgPack/src/SpVectorOp.h"
#include "SparseLinAlgPack/src/MatrixWithOp.h"
#include "LinAlgPack/src/LinAlgOpPack.h"
#include "LinAlgPack/src/VectorOp.h"
#include "LinAlgPack/src/VectorOut.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

bool ReducedSpaceSQPPack::CalcLambdaIndepStd_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{

	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;

	using LinAlgPack::Vt_S;
	using LinAlgPack::norm_inf;

	using SparseLinAlgPack::Vp_StMtV;
	
	using ConstrainedOptimizationPack::min_abs;

	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_StMtV;
	using LinAlgOpPack::Vp_MtV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	Range1D		decomp	= s.equ_decomp();
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Compute: lambda(decomp) = inv(Gc(decomp)'* Y)' * ( - Y' * (Gf + nu) - U' * lambda(undecomp) )
	// where U = Gc(undecomp)' * Y

	// Must resize lambda here explicitly since we will only be updating a region of it.
	// If lambda(undecomp) has already been updated then lambda will have been resized
	// already but lambda(decomp) will not be initialized yet.
	if( !s.lambda().updated_k(0) ) s.lambda().set_k(0).v().resize( algo.nlp().m() );
	
	VectorSlice lambda_decomp = s.lambda().get_k(0).v()(decomp);
	
	// lambda_decomp_tmp1 = - Y' * (Gf + nu)
	if( algo.nlp().has_bounds() ) {
		// _tmp = Gf + nu
		Vector _tmp = s.Gf().get_k(0)();
		VectorSlice _vs_tmp = _tmp;	// only create this VectorSlice once
		Vp_V( &_vs_tmp, s.nu().get_k(0)() );
		// lambda_decomp_tmp1 = - Y' * _tmp
		V_StMtV( &lambda_decomp, -1.0, s.Y().get_k(0), trans, _vs_tmp );
	}
	else {
		// lambda_decomp__tmp1 = - Y' * Gf
		V_StMtV( &lambda_decomp, -1.0, s.Y().get_k(0), trans, s.Gf().get_k(0)() );
	}

	// lambda_decomp_tmp2 = lambda_decomp_tmp1 - U' * lambda(undecomp)
	if( algo.nlp().r() < algo.nlp().m() ) {
		Range1D undecomp = s.equ_undecomp();
		Vp_StMtV( &lambda_decomp, -1.0, s.U().get_k(0), trans, s.lambda().get_k(0).v()(undecomp) );
	}
	// else lambda(decomp)_tmp2 = lambda(decomp)_tmp1

	// lambda_decomp = inv(Gc(decomp)'* Y)' * lambda_decomp_tmp2
	s.decomp_sys().solve_transAtY( lambda_decomp, trans, &lambda_decomp );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nmax(|lambda_k(equ_decomp)(i)|) = " << norm_inf(lambda_decomp)
			<< "\nmin(|lambda_k(equ_decomp)(i)|) = " << min_abs(lambda_decomp)  << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nlambda_k(equ_decomp) = \n" << lambda_decomp;
	}

	return true;
}

void ReducedSpaceSQPPack::CalcLambdaIndepStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculate the Lagrange multipliers for the decomposed constraints\n"
		<< L << "lambda_k(equ_decomp) = - inv(Gc_k(:,equ_decomp)'*Y_k)\n"
		<< L << "                        * (Y_k'*(Gf_k + nu_k) + U_k'*lambda_k(equ_undecomp))\n";
}
