// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianExactStd_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <sstream>
#include <typeinfo>
#include <iomanip>

#include "ReducedSpaceSQPPack/include/std/ReducedHessianExactStd_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "SparseLinAlgPack/include/MatrixSymDenseInitialize.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "NLPInterfacePack/include/NLPSecondOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/GenMatrixAsTriSym.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

bool ReducedHessianExactStd_Step::do_step(
	  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	, poss_type assoc_step_poss)
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgPack::nonconst_sym;
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
		nlp.get_lambda_init( &s.lambda().set_k(-1).v() );
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "||lambda_km1||inf = " << s.lambda().get_k(-1).norm_inf() << std::endl;
		}
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
			out << "lambda_km1 = \n" << s.lambda().get_k(-1)();
		}
	}

	nlp.set_HL(	&s.HL().set_k(0) );
	nlp.calc_HL( s.x().get_k(0)(), s.lambda().get_k(-1)(), false );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
		s.HL().get_k(0).output( out << "\nHL_k = \n" );
	}

	// If rHL has already been updated for this iteration then just leave it.
	if( !s.rHL().updated_k(0) ) {

		const MatrixSymWithOp
			*HL_sym_op = dynamic_cast<const MatrixSymWithOp*>(&s.HL().get_k(0));
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
		GenMatrix rHL_sym_store(nind,nind);
		sym_gms rHL_sym(rHL_sym_store(),BLAS_Cpp::lower);
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
