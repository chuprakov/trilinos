// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianExactStd_Step.cpp
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

#include <sstream>
#include <typeinfo>
#include <iomanip>

#include "ReducedSpaceSQPPack/src/std/ReducedHessianExactStd_Step.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "SparseLinAlgPack/src/MatrixSymDenseInitialize.hpp"
#include "GeneralIterationPack/src/print_algorithm_step.hpp"
#include "ConstrainedOptimizationPack/src/VectorWithNorms.h"
#include "NLPInterfacePack/src/NLPSecondOrderInfo.hpp"
#include "SparseLinAlgPack/src/MatrixSymWithOp.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "DenseLinAlgPack/src/DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack/src/DMatrixOut.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DVectorOp.hpp"
#include "DenseLinAlgPack/src/DVectorOut.hpp"
#include "Midynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

bool ReducedHessianExactStd_Step::do_step(
	  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	, poss_type assoc_step_poss)
{
	using DynamicCastHelperPack::dyn_cast;
	using DenseLinAlgPack::nonconst_sym;
	using SparseLinAlgPack::Mp_StMtMtM;
	typedef SparseLinAlgPack::MatrixSymDenseInitialize	MatrixSymDenseInitialize;
	typedef SparseLinAlgPack::MatrixSymWithOp			MatrixSymWithOp;
	using ConstrainedOptimizationPack::NLPSecondOrderInfo;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLPSecondOrderInfo
#ifdef _WINDOWS
				&nlp	= dynamic_cast<NLPSecondOrderInfo&>(algo.nlp());
#else
				&nlp	= dyn_cast<NLPSecondOrderInfo>(algo.nlp());
#endif
	MatrixSymWithOp
		*HL_sym_op = dynamic_cast<MatrixSymWithOp*>(&s.HL().get_k(0));

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// problem size
	size_type	n		= nlp.n(),
				r		= nlp.r(),
				nind	= n - r;

	// Compute HL first (You may want to move this into its own step later?)

	if( !s.lambda().updated_k(-1) ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "Initializing lambda_km1 = nlp.get_lambda_init ... \n";
		}
		nlp.get_init_lagrange_mult( &s.lambda().set_k(-1).v(), NULL );
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "||lambda_km1||inf = " << s.lambda().get_k(-1).norm_inf() << std::endl;
		}
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
			out << "lambda_km1 = \n" << s.lambda().get_k(-1)();
		}
	}

	nlp.set_HL(	HL_sym_op );
	nlp.calc_HL( s.x().get_k(0)(), s.lambda().get_k(-1)(), false );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		s.HL().get_k(0).output( out << "\nHL_k = \n" );
	}

	// If rHL has already been updated for this iteration then just leave it.
	if( !s.rHL().updated_k(0) ) {

		if( !HL_sym_op ) {
			std::ostringstream omsg;
			omsg
				<< "ReducedHessianExactStd_Step::do_step(...) : Error, "
				<< "The matrix HL with the concrete type "
				<< typeid(s.HL().get_k(0)).name() << " does not support the "
				<< "MatrixSymWithOp iterface";
			throw std::logic_error( omsg.str() );
		}		

		MatrixSymDenseInitialize
			*rHL_sym_init = dynamic_cast<MatrixSymDenseInitialize*>(&s.rHL().set_k(0));
		if( !rHL_sym_init ) {
			std::ostringstream omsg;
			omsg
				<< "ReducedHessianExactStd_Step::do_step(...) : Error, "
				<< "The matrix rHL with the concrete type "
				<< typeid(s.rHL().get_k(0)).name() << " does not support the "
				<< "MatrixSymDenseInitialize iterface";
			throw std::logic_error( omsg.str() );
		}		

		// Compute the dense reduced Hessian
		DMatrix rHL_sym_store(nind,nind);
		DMatrixSliceSym rHL_sym(rHL_sym_store(),BLAS_Cpp::lower);
		Mp_StMtMtM( &rHL_sym, 1.0, MatrixSymWithOp::DUMMY_ARG, *HL_sym_op
					, s.Z().get_k(0), BLAS_Cpp::no_trans, 0.0 );

		if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
			out << "\nLower triangular partion of dense reduced Hessian (ignore nonzeros above diagonal):\nrHL_dense = \n" << rHL_sym_store(); 
		}
	
		// Set the reduced Hessain
		rHL_sym_init->initialize( rHL_sym );

		if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
			s.rHL().get_k(0).output( out << "\nrHL_k = \n" );
		}
	}

	return true;
}

void ReducedHessianExactStd_Step::print_step(
	  const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	, poss_type assoc_step_poss, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculate the exact reduced Hessian of the Lagrangian\n"
		<< L << "if lambda_km1 is not updated then\n"
		<< L << "    lambda_km1 = nlp.get_lambda_init\n"
		<< L << "end\n"
		<< L << "HL_k = HL(x_k,lambda_km1) <: R^(n+m) -> R^(n x n)\n"
		<< L << "if rHL_k is not updated then\n"
		<< L << "    rHL_dense = Z_k' * HL_k * Z_k  (MatrixSymWithOp interface for HL_k)\n"
		<< L << "    rHL_k = rHL_dense (MatrixSymDenseInitialize interface for rHL_k)\n"
		<< L << "end\n";
}

}	// end namespace ReducedSpaceSQPPack 
