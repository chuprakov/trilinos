// ////////////////////////////////////////////////////////////////////////////
// NullSpaceStepIP_Step.cpp
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
#include <iostream>

#include "ReducedSpaceSQPPack/include/std/NullSpaceStepIP_Step.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_Step.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

bool NullSpaceStepIP_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
	{
	using BLAS_Cpp::no_trans;
	using DynamicCastHelperPack::dyn_cast;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using AbstractLinAlgPack::Vt_S;
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::V_InvMtV;
	using AbstractLinAlgPack::Mp_StM;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_MtV;
 	using LinAlgOpPack::M_StM;
	using LinAlgOpPack::assign;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	ipState	    &s      = dyn_cast<ipState>(_algo.state());

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Compute qp_grad which is an approximation to rGf + Z'*(mu*(invXu*e-invXl*e) + no_cross_term
	// minimize round off error by calc'ing Z'*(Gf + mu*(invXu*e-invXl*e))

	// qp_grad_k = Z'*(Gf + mu*(invXu*e-invXl*e))
	const MatrixSymDiagonalStd& invXu = s.invXu().get_k(0);
	const MatrixSymDiagonalStd& invXl = s.invXl().get_k(0);
	const value_type& mu = s.barrier_parameter().get_k(0);

	MemMngPack::ref_count_ptr<VectorWithOpMutable> rhs = s.Gf().get_k(0).clone();
	rhs->axpy(mu, invXu.diag());
	rhs->axpy(-1.0*mu, invXl.diag());
	
	out << "mu(invXu-invXl)\n";
	rhs->output(out);

	VectorWithOpMutable &qp_grad_k = s.qp_grad().set_k(0);
	V_MtV(&qp_grad_k, s.Z().get_k(0), BLAS_Cpp::trans, *rhs);
	
	out << "rhs (no w)=\n";
	qp_grad_k.output(out);

	// error check for cross term
	value_type& zeta = s.zeta().set_k(0);
	const VectorWithOp& w_sigma = s.w_sigma().get_k(0);
	
	// need code to calculate damping parameter
	zeta = 1.0;

	using LinAlgOpPack::Vp_StV;
	Vp_StV(&qp_grad_k, zeta, w_sigma);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "rhs=\n";
		qp_grad_k.output(out);
		}

	// build the "Hessian" term B = rHL + rHB
	// should this be MatrixSymWithOpNonsingular
	const MatrixSymWithOpNonsingular  &rHL_k = dyn_cast<MatrixSymWithOpNonsingular>(s.rHL().get_k(0));
	const MatrixSymWithOp  &rHB_k = dyn_cast<MatrixSymWithOpNonsingular>(s.rHB().get_k(0));
	const MatrixWithOp& Z_k = s.Z().get_k(0);

	MatrixSymWithOpNonsingular& B_k = dyn_cast<MatrixSymWithOpNonsingular>(s.B().set_k(0));
	if (B_k.cols() != Z_k.cols())
		{
		// Initialize space in rHB
		dyn_cast<MatrixSymInitDiagonal>(B_k).init_identity(Z_k.space_rows(), 0.0);
		}

	//	M_StM(&B_k, 1.0, rHL_k, no_trans);
	assign(&B_k, rHL_k, BLAS_Cpp::no_trans);
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "rHL_k=\n";
		B_k.output(out);
		}

	Mp_StM(&B_k, 1.0, rHB_k, BLAS_Cpp::no_trans);
	
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "B=\n";
		B_k.output(out);
		}
	

	// Solve the system pz = - inv(rHL) * qp_grad
	VectorWithOpMutable               &pz_k  = s.pz().set_k(0);
	V_InvMtV( &pz_k, B_k, no_trans, qp_grad_k );
	Vt_S( &pz_k, -1.0 );

	// Zpz = Z * pz
	V_MtV( &s.Zpz().set_k(0), s.Z().get_k(0), no_trans, pz_k );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||pz||inf   = " << s.pz().get_k(0).norm_inf()
			<< "\nsum(Zpz)    = " << AbstractLinAlgPack::sum(s.Zpz().get_k(0))  << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\npz_k = \n" << s.pz().get_k(0);
		out << "\nnu_k = \n" << s.nu().get_k(0);
		out << "\nZpz_k = \n" << s.Zpz().get_k(0);
		out << std::endl;
	}

	if(algo.algo_cntr().check_results()) {
		assert_print_nan_inf(s.pz().get_k(0),  "pz_k",true,&out);
		assert_print_nan_inf(s.Zpz().get_k(0), "Zpz_k",true,&out);
	}

	return true;
	}

void NullSpaceStepIP_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
	{
	out
		<< L << "*** Calculate the null space step by solving an unconstrainted QP\n"
		<< L << "TODO: Correct this documentation\n";
	}

} // end namespace ReducedSpaceSQPPack
