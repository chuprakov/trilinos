// ////////////////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateStd_Step.cpp
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

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateStd_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixSymInitDiagonal.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
}

ReducedSpaceSQPPack::ReducedHessianSecantUpdateStd_Step::ReducedHessianSecantUpdateStd_Step(
	const secant_update_ptr_t&   secant_update
	)
	:secant_update_(secant_update)
	 ,iter_k_rHL_init_ident_(-1.0) // not a valid iteration?
{}

bool ReducedSpaceSQPPack::ReducedHessianSecantUpdateStd_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	using DynamicCastHelperPack::dyn_cast;
	using AbstractLinAlgPack::dot;
	using AbstractLinAlgPack::Vt_S;
	using AbstractLinAlgPack::Vp_StV;
	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_VpV;
	using LinAlgOpPack::V_VmV;
	using LinAlgOpPack::V_StV;
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

	bool return_val = true;
	
	// Get iteration quantities
	IterQuantityAccess<index_type>
		&num_basis_iq = s.num_basis();
	IterQuantityAccess<VectorWithOpMutable>
		&py_iq  = s.py(),
		&pz_iq  = s.pz(),
		&rGf_iq = s.rGf(),
		&w_iq   = s.w();
	IterQuantityAccess<MatrixWithOp>
		&Z_iq = s.Z();
	IterQuantityAccess<MatrixSymWithOp>
		&rHL_iq = s.rHL();

	// problem size
	const size_type
		n    = algo.nlp().n(),
		r    = py_iq.get_k(py_iq.last_updated()).dim(),
		nind = n - r;

	// See if a new basis has been selected
	bool new_basis = false;
	{
		const int last_updated_k = num_basis_iq.last_updated();
		if( last_updated_k != IterQuantity::NONE_UPDATED ) {
			const index_type num_basis_last = num_basis_iq.get_k(last_updated_k);
			if( num_basis_ == NO_BASIS_UPDATED_YET )
				num_basis_ = num_basis_last;
			else if( num_basis_ != num_basis_last )
				new_basis = true;
			num_basis_ = num_basis_last;
		}
	}

	// If rHL has already been updated for this iteration then just leave it.
	if( !rHL_iq.updated_k(0) ) {

		// If a new basis has been selected, reinitialize
		if( new_basis ) {

			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nBasis changed.  Reinitializing rHL_k = eye(n-r) ...\n";
			}
			dyn_cast<MatrixSymInitDiagonal>(rHL_iq.set_k(0)).init_identity(
				Z_iq.get_k(0).space_rows()
				);
			if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES )
				out << "\nrHL_k = \n" << rHL_iq.get_k(0);
			quasi_newton_stats_(s).set_k(0).set_updated_stats(
				QuasiNewtonStats::REINITIALIZED );
			iter_k_rHL_init_ident_ = s.k();	// remember what iteration this was

		}
		else {
			
			// Determine if rHL has been initialized and if we
			// can perform the update.  To perform the BFGS update
			// rHL_km1 and rGf_km1 must have been computed.
			if( rHL_iq.updated_k(-1) && rGf_iq.updated_k(-1) ) {

				// /////////////////////////////////////////////////////
				// Perform the Secant update

				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
				{
					out
						<< "\nPerforming Secant update ...\n";
				}

				const VectorWithOp
					&rGf_k   = rGf_iq.get_k(0),
					&rGf_km1 = rGf_iq.get_k(-1),
					&pz_km1  = pz_iq.get_k(-1);
				const value_type
					alpha_km1 = s.alpha().get_k(-1);
				VectorSpace::vec_mut_ptr_t
					y_bfgs = rGf_k.space().create_member(),
					s_bfgs = pz_km1.space().create_member();

				// /////////////////////////////////////////////////////
				// y_bfgs = rGf_k - rGf_km1 - alpha_km1 * w_km1
			
				// y_bfgs = rGf_k - rGf_km1 
				V_VmV( y_bfgs.get(), rGf_k, rGf_km1 );	

				if( w_iq.updated_k(-1) )
					// y_bfgs += - alpha_km1 * w_km1
					Vp_StV( y_bfgs.get(), - alpha_km1, w_iq.get_k(-1) );

				// /////////////////////////////////////////////////////
				// s_bfgs = alpha_km1 * pz_km1
				V_StV( s_bfgs.get(), alpha_km1, pz_iq.get_k(-1) );

				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
					out << "\n||y_bfgs||inf = " << y_bfgs->norm_inf() << std::endl;
					out << "\n||s_bfgs||inf = " << s_bfgs->norm_inf() << std::endl;
				}

				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
					out << "\ny_bfgs =\n" << *y_bfgs;
					out << "\ns_bfgs =\n" << *s_bfgs;
				}

				// Update from last
				MatrixSymWithOp
					&rHL_k   = rHL_iq.set_k(0,-1);

				// Perform the secant update
				if(!secant_update().perform_update(
					   s_bfgs.get(), y_bfgs.get()
					   ,iter_k_rHL_init_ident_ == s.k() - 1
					   ,out, olevel, &algo, &s, &rHL_k
					   ))
				{
					return_val = false; // redirect control of algorithm!
				}

			}
			else {
				// We do not have the info to perform the update

				int k_last_offset = rHL_iq.last_updated();
				bool set_current = false;
				if( k_last_offset != IterQuantity::NONE_UPDATED && k_last_offset < 0 ) {
					const MatrixSymWithOp &rHL_k_last = rHL_iq.get_k(k_last_offset);
					const size_type nind_last = rHL_k_last.rows();
					if( nind_last != nind) {
						if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
							out
								<< "No new basis was selected.\n"
								<< "The previous matrix rHL_k(" << k_last_offset << ") was found but its dimmension\n"
								<< "rHL_k(" << k_last_offset << ").rows() = " << nind_last << " != n-r = " << nind << std::endl;
						}
					}
					else {
						if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
							out
								<< "No new basis was selected so using previously updated...\n "
								<< "rHL_k = rHL_k(" << k_last_offset << ")\n";
						}
						rHL_iq.set_k(0) = rHL_k_last;
						quasi_newton_stats_(s).set_k(0).set_updated_stats(
							QuasiNewtonStats::SKIPED );
					
						if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
							rHL_iq.get_k(0).output( out << "\nrHL_k = \n" );
						}
						set_current = true;
					}
				}
				if( !set_current ) {
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out
							<< "\nInitializing rHL = eye(n-r) "
							<< "(k = " << algo.state().k() << ")...\n";
					}

					// Now I will assume that since I can't perform the BFGS update and rHL has
					// not been set for this iteration yet, that it is up to me to initialize rHL_k = 0
					dyn_cast<MatrixSymInitDiagonal>(rHL_iq.set_k(0)).init_identity(
						Z_iq.get_k(0).space_rows() );
					iter_k_rHL_init_ident_ = s.k();	// remember what iteration this was
					quasi_newton_stats_(s).set_k(0).set_updated_stats(
						QuasiNewtonStats::REINITIALIZED );
				
				}
			}

		}
		
		// Print rHL_k
		
		MatrixWithOp::EMatNormType mat_nrm_inf = MatrixWithOp::MAT_NORM_INF;

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\n||rHL_k||inf    = " << rHL_iq.get_k(0).calc_norm(mat_nrm_inf).value << std::endl;
			if(algo.algo_cntr().calc_conditioning()) {
				const MatrixSymWithOpNonsingular
					*rHL_ns_k = dynamic_cast<const MatrixSymWithOpNonsingular*>(&rHL_iq.get_k(0));
				if(rHL_ns_k)
					out << "\ncond_inf(rHL_k) = " << rHL_ns_k->calc_cond_num(mat_nrm_inf).value << std::endl;
			}
		}

		if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
			out << "\nrHL_k = \n" << rHL_iq.get_k(0);
		}
		
	}
	else {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "\nThe matrix rHL_k has already been updated so leave it\n";
		}
	}
	
	return return_val;
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
		<< L << "        rHL_k = rHL_km1\n"
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
