// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateStd_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <ostream>

#include "../../include/std/ReducedHessianSecantUpdateStd_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "../../include/ReducedSpaceSQPPackExceptions.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/test/TestMatrixSymSecantUpdate.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/MatrixSymInitDiagonal.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

ReducedSpaceSQPPack::ReducedHessianSecantUpdateStd_Step::ReducedHessianSecantUpdateStd_Step(
		  const secant_update_ptr_t&   secant_update
		)
	:
	secant_update_(secant_update)
	,iter_k_rHL_init_ident_(-1.0) // not a valid iteration?
{}

bool ReducedSpaceSQPPack::ReducedHessianSecantUpdateStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgPack::norm_2;
	using LinAlgPack::norm_inf;
	using LinAlgPack::dot;
	using LinAlgPack::Vt_S;
	using LinAlgPack::V_VpV;
	using LinAlgPack::V_VmV;
	using LinAlgPack::Vp_StV;
	using LinAlgPack::V_StV;

	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_MtV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// problem size
	size_type	n		= algo.nlp().n(),
				r		= algo.nlp().r(),
				nind	= n - r;

	// Initialize basis counter
	if( num_basis_ == NO_BASIS_UPDATED_YET ) num_basis_ = s.num_basis();

	// See if a new basis has been selected
	bool new_basis;
	if( new_basis = ( num_basis_ != s.num_basis() ) )
		num_basis_ = s.num_basis();

	// If rHL has already been updated for this iteration then just leave it.
	if( !s.rHL().updated_k(0) ) {

		// If a new basis has been selected, reinitialize
		if( new_basis ) {

			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nBasis changed.  Reinitializing rHL_k = eye(n-r)\n";
			}
#ifdef _WINDOWS
			dynamic_cast<MatrixSymInitDiagonal&>(s.rHL().set_k(0)).init_identity(nind);
#else
			dyn_cast<MatrixSymInitDiagonal>(s.rHL().set_k(0)).init_identity(nind);
#endif
			quasi_newton_stats_(s).set_k(0).set_updated_stats(
				QuasiNewtonStats::REINITIALIZED );
			iter_k_rHL_init_ident_ = s.k();	// remember what iteration this was
			return true;
		}

		// Determine if rHL has been initialized and if we
		// can perform the update.  To perform the BFGS update
		// rHL_km1 and rGf_km1 must have been computed.
		if( s.rHL().updated_k(-1) && s.rGf().updated_k(-1) ) {

			// /////////////////////////////////////////////////////
			// Perform the Secant update

			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
			{
				out
					<< "\nPerforming Secant update ...\n";
			}

			Vector y_bfgs, s_bfgs;

			// /////////////////////////////////////////////////////
			// y_bfgs = rGf_k - rGf_km1 - alpha_km1 * w_km1
			
			// y_bfgs = rGf_k - rGf_km1 
			V_VmV( &y_bfgs, s.rGf().get_k(0)(), s.rGf().get_k(-1)() );	

			if( s.w().updated_k(-1) )
				// y_bfgs += - alpha_km1 * w_km1
				Vp_StV( &y_bfgs(), - s.alpha().get_k(-1), s.w().get_k(-1)() );

			// /////////////////////////////////////////////////////
			// s_bfgs = alpha_km1 * pz_km1
			V_StV( &s_bfgs, s.alpha().get_k(-1), s.pz().get_k(-1)() );

			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << "\n||y_bfgs||inf = " << norm_inf(y_bfgs()) << std::endl;
				out << "\n||s_bfgs||inf = " << norm_inf(s_bfgs()) << std::endl;
			}

			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
				out << "\ny_bfgs =\n" << y_bfgs;
				out << "\ns_bfgs =\n" << s_bfgs;
			}

			// Update from last
			MatrixWithOp
				&rHL_km1 = s.rHL().get_k(-1),
				&rHL_k = s.rHL().set_k(0) = rHL_km1;

			// Perform the secant update
			if(!secant_update().perform_update(
				&s_bfgs(), &y_bfgs(), iter_k_rHL_init_ident_ == s.k() - 1
				,out, olevel, &algo, &s, &rHL_k
				))
			{
				return false; // redirect control of algorithm!
			}

		}
		else {
			// We do not have the info to perform the update

			int k_last_offset = s.rHL().last_updated();

			if( k_last_offset == IterQuantity::NONE_UPDATED && k_last_offset < 0 ) {

				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
					algo.track().journal_out()
						<< "\nInitializing rHL = eye(n-r) "
						<< "(k = " << algo.state().k() << ")...\n";
				}

				// Now I will assume that since I can't perform the BFGS update and rHL has
				// not been set for this iteration yet, that it is up to me to initialize rHL_k = 0
#ifdef _WINDOWS
				dynamic_cast<MatrixSymInitDiagonal&>(s.rHL().set_k(0)).init_identity(nind);
#else
				dyn_cast<MatrixSymInitDiagonal>(s.rHL().set_k(0)).init_identity(nind);
#endif
				iter_k_rHL_init_ident_ = s.k();	// remember what iteration this was
				quasi_newton_stats_(s).set_k(0).set_updated_stats(
					QuasiNewtonStats::REINITIALIZED );

				if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
					s.rHL().get_k(0).output( out << "\nrHL_k = \n" );
				}

			}
			else {

				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
					algo.track().journal_out()
						<< "No new basis was selected so using previously updated...\n "
						<< "rHL_k = rHL_k(" << k_last_offset << ")\n";
				}
				const MatrixWithOp &rHL_k_last = s.rHL().get_k(k_last_offset);
				s.rHL().set_k(0) = rHL_k_last;
				quasi_newton_stats_(s).set_k(0).set_updated_stats(
					QuasiNewtonStats::SKIPED );

				if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
					s.rHL().get_k(0).output( out << "\nrHL_k = \n" );
				}

			}
		}
	}
	else {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "The matrix rHL_k has already been updated.\n";
		}
	}

	return true;
}

void ReducedSpaceSQPPack::ReducedHessianSecantUpdateStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculate the reduced hessian of the Lagrangian rHL = Z' * HL * Z\n"
		<< L << "default:  num_basis_remembered = NO_BASIS_UPDATED_YET\n"
		<< L << "          iter_k_rHL_init_ident = -1\n"
		<< L << "if num_basis_remembered = NO_BASIS_UPDATED_YET then\n"
		<< L << "    num_basis_remembered = num_basis\n"
		<< L << "end\n"
		<< L << "if num_basis_remembered != num_basis then\n"
		<< L << "    num_basis_remembered = num_basis\n"
		<< L << "    new_basis = true\n"
		<< L << "end\n"
		<< L << "if rHL_k is not updated then\n"
		<< L << "    if new_basis == true then\n"
		<< L << "        *** Transition rHL to the new basis by just starting over.\n"
		<< L << "        rHL_k = eye(n-r) *** must support MatrixSymInitDiagonal interface\n"
		<< L << "        iter_k_rHL_init_ident = k\n"
		<< L << "        goto next step\n"
		<< L << "    end\n"
		<< L << "    if rHL_km1 and rGf_km1 are updated then\n"
		<< L << "        *** We should have the information to perform a BFGS update\n"
		<< L << "        y = rGf_k - rGf_km1\n"
		<< L << "        s = alpha_km1 * pz_km1\n"
		<< L << "        if k - 1 == iter_k_rHL_init_ident then\n"
		<< L << "            first_update = true\n"
		<< L << "        else\n"
		<< L << "            first_update = false\n"
		<< L << "        end\n"
		<< L << "        begin secant update\n"
		<< L << "        (" << typeid(secant_update()).name() << ")\n"
		;
	secant_update().print_step( out, L+"            " );
	out
		<< L << "        end secant update\n"
		<< L << "    else\n"
		<< L << "       *** We have no information for which to preform a BFGS update.\n"
		<< L << "       k_last_offset = last iteration rHL was updated for\n"
		<< L << "       if k_last_offset does not exist then\n"
		<< L << "            *** We are left with no choise but to initialize rHL\n"
		<< L << "            rHL_k = eye(n-r) *** must support MatrixSymInitDiagonal interface\n"
		<< L << "            iter_k_rHL_init_ident = k\n"
		<< L << "        else\n"
		<< L << "            *** No new basis has been selected so we may as well\n"
		<< L << "            *** just use the last rHL that was updated\n"
		<< L << "            rHL_k = rHL_k(k_last_offset)\n"
		<< L << "        end\n"
		<< L << "    end\n"
		<< L << "end\n";
}
