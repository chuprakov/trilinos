// ////////////////////////////////////////////////////////////////////////////
// CalcLambdaIndepStd_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "../../include/std/CalcLambdaIndepStd_AddedStep.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/ComputeMinMult.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"

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
	Range1D		indep	= s.con_indep();
	
	EIterationInfoOutput olevel = s.iteration_info_output();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Compute: lambda(indep) = inv(Gc(indep)'* Y)' * ( - Y' * (Gf + nu) - U' * lambda(dep) )
	// where U = Gc(dep)' * Y

	// Must resize lambda here explicitly since we will only be updating a region of it.
	// If lambda(dep) has already been updated then lambda will have been resized
	// already but lambda(indep) will not be initialized yet.
	if( !s.lambda().updated_k(0) ) s.lambda().set_k(0).v().resize( algo.nlp().m() );
	
	VectorSlice lambda_indep = s.lambda().get_k(0).v()(indep);
	
	// lambda_indep_tmp1 = - Y' * (Gf + nu)
	if( algo.nlp().has_bounds() ) {
		// _tmp = Gf + nu
		Vector _tmp = s.Gf().get_k(0)();
		VectorSlice _vs_tmp = _tmp;	// only create this VectorSlice once
		Vp_V( &_vs_tmp, s.nu().get_k(0)() );
		// lambda_indep_tmp1 = - Y' * _tmp
		V_StMtV( &lambda_indep, -1.0, s.Y().get_k(0), trans, _vs_tmp );
	}
	else {
		// lambda_indep__tmp1 = - Y' * Gf
		V_StMtV( &lambda_indep, -1.0, s.Y().get_k(0), trans, s.Gf().get_k(0)() );
	}

	// lambda_indep_tmp2 = lambda_indep_tmp1 - U' * lambda(dep)
	if( algo.nlp().r() < algo.nlp().m() ) {
		Range1D dep = s.con_dep();
		Vp_StMtV( &lambda_indep, -1.0, s.U().get_k(0), trans, s.lambda().get_k(0).v()(dep) );
	}
	// else lambda(indep)_tmp2 = lambda(indep)_tmp1

	// lambda_indep = inv(Gc(indep)'* Y)' * lambda_indep_tmp2
	s.decomp_sys().solve_transAtY( lambda_indep, trans, &lambda_indep );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nmax(|lambda_k(con_indep)(i)|) = " << norm_inf(lambda_indep)
			<< "\nmin(|lambda_k(con_indep)(i)|) = " << min_abs(lambda_indep)  << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\nlambda_k(con_indep) = \n" << lambda_indep;
	}

	return true;
}

void ReducedSpaceSQPPack::CalcLambdaIndepStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculate the Lagrange multipliers for the independent constraints\n"
		<< L << "lambda_k(con_indep) = - inv(Gc_k(:,con_indep)'*Y_k)\n"
		<< L << "                        * (Y_k'*(Gf_k + nu_k) + U_k'*lambda_k(con_dep))\n";
}
