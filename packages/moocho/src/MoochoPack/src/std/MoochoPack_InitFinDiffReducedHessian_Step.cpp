// ////////////////////////////////////////////////////////////////////////////
// InitFinDiffReducedHessian_Step.cpp
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

#include "ReducedSpaceSQPPack/include/std/InitFinDiffReducedHessian_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "NLPInterfacePack/include/NLPObjGradient.h"
#include "SparseLinAlgPack/include/MatrixSymInitDiagonal.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/max_near_feas_step.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

ReducedSpaceSQPPack::InitFinDiffReducedHessian_Step::InitFinDiffReducedHessian_Step(
		  EInitializationMethod		initialization_method
		, value_type				max_cond
		, value_type				min_diag
		, value_type				step_scale			)
	: num_basis_(NO_BASIS_UPDATED_YET)
		, initialization_method_(initialization_method)
		, max_cond_(max_cond), min_diag_(min_diag), step_scale_(step_scale)
{}

bool ReducedSpaceSQPPack::InitFinDiffReducedHessian_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgPack::norm_inf;
	using LinAlgPack::Vt_S;
	using LinAlgPack::Vp_StV;
	using LinAlgPack::norm_inf;
	using LinAlgOpPack::V_MtV;
	using SparseLinAlgPack::max_near_feas_step;
	using NLPInterfacePack::NLPObjGradient;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
#ifdef _WINDOWS
	NLPObjGradient
		&nlp = dynamic_cast<NLPObjGradient&>(algo.nlp());
#else
	NLPObjGradient
		&nlp = dyn_cast<NLPObjGradient>(algo.nlp());
#endif
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Initialize basis counter
	if( num_basis_ == NO_BASIS_UPDATED_YET ) num_basis_ = s.num_basis();

	// See if a new basis has been selected
	bool new_basis;
	if( new_basis = ( num_basis_ != s.num_basis() ) )
		num_basis_ = s.num_basis();

	int k_last_offset = s.rHL().last_updated();

	// If the basis has changed or there is no previous matrix to use
	// then reinitialize.

	if( new_basis || k_last_offset == IterQuantity::NONE_UPDATED ) {

		// Compute a finite difference along the null space of the
		// constraints

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\nReinitializing the reduced Hessain using a finite difference\n";
		}

#ifdef _WINDOWS
		MatrixSymInitDiagonal
			&rHL_diag = dynamic_cast<MatrixSymInitDiagonal&>(s.rHL().set_k(0));
#else
		MatrixSymInitDiagonal
			&rHL_diag = dyn_cast<MatrixSymInitDiagonal>(s.rHL().set_k(0));
#endif

		// problem size
		size_type	n		= algo.nlp().n(),
					r		= algo.nlp().r(),
					nind	= n - r;

		// one vector
		Vector e(nind);
		e = 1.0;

		// Ze
		Vector Ze;
		V_MtV( &Ze, s.Z().get_k(0), BLAS_Cpp::no_trans, e() );
		
		// This does not have to be an accurate finite difference so lets just
		// take step_scale/||Ze|| as the step size unless it is out of bounds
		// If we assume that our variables are scaled near
		// one (step_scale == 1?) then this will give us an appreciable step.  Beside we
		// should not be near the solution so the reduced gradient should not
		// be near zero.

		const value_type nrm_Ze = norm_inf( Ze() );
		value_type u = step_scale() / nrm_Ze;

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\n||Ze||inf = " << nrm_Ze << std::endl;
		}

		if( (int)olevel >= (int)PRINT_VECTORS ) {
			out << "\nZe =\n" << Ze();
		}

		if( algo.nlp().has_bounds() ) {

			// Find the maximum step u
			// we can take along x_fd = x_k + u*Ze
			// that don't violate variable bounds by too much.
			// If a positive step can't be found then this may be a negative step.
			
			std::pair<value_type,value_type>
				u_steps = max_near_feas_step( s.x().get_k(0)(), Ze()
					, algo.nlp().xl(), algo.nlp().xu()
					, algo.algo_cntr().max_var_bounds_viol() );
			
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nMaximum steps ( possitive, negative ) in bounds u = ("
					<< u_steps.first << "," << u_steps.second << ")\n";
			}

			if( u_steps.first < u )
				u = u_steps.first;
			if( ::fabs(u_steps.second) < u )
				u = u_steps.second;
		}

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\nFinite difference step length u = " << u << std::endl;
		}

		// Take the finite difference from x_k to x_fd = x_k + u * Ze
		//
		// rGf_fd = ( Z_k'*Gf(x_k + u*Ze) - rGf_k ) / u
		//

		Vector x_fd = s.x().get_k(0)();
		Vp_StV( &x_fd(), u, Ze() );

		// Gf_fd = Gf(x_fd)
		Vector Gf_fd;
		nlp.set_Gf(	&Gf_fd );
		nlp.set_multi_calc( false );
		nlp.calc_Gf( x_fd );

		if( (int)olevel >= (int)PRINT_VECTORS ) {
			out << "\nGf_fd =\n" << Gf_fd();
		}

		// rGf_fd = Z'*Gf_fd
		Vector rGf_fd;
		V_MtV( &rGf_fd, s.Z().get_k(0), BLAS_Cpp::trans, Gf_fd() );

		// rGf_fd = rGf_fd - rGf_k
		Vp_StV( &rGf_fd(), -1.0, s.rGf().get_k(0)() );

		// rGf_fd = rGf_fd / u
		Vt_S( &rGf_fd(), 1.0 / u );

		const value_type
			nrm_rGf_fd = norm_inf( rGf_fd() );

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\n||(rGf_fd - rGf_k)/u||inf = " << nrm_rGf_fd << std::endl;
		}
		if( (int)olevel >= (int)PRINT_VECTORS ) {
			out << "\n(rGf_fd - rGf_k)/u =\n" << rGf_fd();
		}

		if( nrm_rGf_fd <= min_diag() ) {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\n||(rGf_fd - rGf_k)/u||inf = " << nrm_rGf_fd
						<< " < min_diag = " << min_diag() << std::endl
					<< "\nScale by min_diag ... \n";
			}
			rHL_diag.init_identity(nind,min_diag());
		}
		else {
			switch( initialization_method() ) {
				case SCALE_IDENTITY: {
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out << "\nScale the identity matrix by ||(rGf_fd - rGf_k)/u||inf ... \n";
					}
					rHL_diag.init_identity(nind,nrm_rGf_fd);
					break;
				}
				case SCALE_DIAGONAL:
				case SCALE_DIAGONAL_ABS:
				{
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out << "\nScale diagonal by modified finite difference ... \n";
					}
					// In order to keep the initial reduced Hessian well conditioned we
					// will not let any diagonal element drop below
					// ||rGf_fd||inf / max_cond
				
					const value_type
						min_ele = std::_MAX( nrm_rGf_fd / max_cond(), min_diag() );

					if( initialization_method() == SCALE_DIAGONAL ) {
						for( Vector::iterator itr = rGf_fd.begin(); itr != rGf_fd.end(); ++itr )
							*itr = std::_MAX( *itr, min_ele );
					}
					else {
						for( Vector::iterator itr = rGf_fd.begin(); itr != rGf_fd.end(); ++itr )
							*itr = std::_MAX( ::fabs(*itr), min_ele );
					}
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out << "\n||diag||inf = " << norm_inf(rGf_fd()) << std::endl;
					}
					if( (int)olevel >= (int)PRINT_VECTORS ) {
						out << "\ndiag =\n" << rGf_fd();
					}
					rHL_diag.init_diagonal(rGf_fd());
					break;
				}
				default:
					assert(0);	// only local programming error?
			}
		}

		quasi_newton_stats_(s).set_k(0).set_updated_stats(
			QuasiNewtonStats::REINITIALIZED );

		if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
			s.rHL().get_k(0).output( out << "\nrHL_k = \n" );
		}

	}

	return true;
}

void ReducedSpaceSQPPack::InitFinDiffReducedHessian_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Initialize the reduced Hessian using a single finite difference.\n"
		<< L << "*** Where the nlp must support the NLPObjGradient interface and\n"
		<< L << "*** rHL_k must support the MatrixSymInitDiagonal interface or exceptions\n"
		<< L << "*** will be thrown.\n"
		<< L << "default: num_basis_remembered = NO_BASIS_UPDATED_YET\n"
		<< L << "         initialization_method = SCALE_DIAGONAL\n"
		<< L << "         max_cond              = " << max_cond() << std::endl
		<< L << "         min_diag              = " << min_diag() << std::endl
		<< L << "         step_scale            = " << step_scale() << std::endl
		<< L << "if num_basis_remembered = NO_BASIS_UPDATED_YET then\n"
		<< L << "    num_basis_remembered = num_basis\n"
		<< L << "end\n"
		<< L << "if num_basis_remembered != num_basis then\n"
		<< L << "    num_basis_remembered = num_basis\n"
		<< L << "    new_basis = true\n"
		<< L << "end\n"
		<< L << "if new_basis == true or no past rHL as been updated then\n"
		<< L << "    *** Reinitialize the reduced Hessian using finite differences\n"
		<< L << "    Ze = Z * e\n"
		<< L << "    u = step_scale / norm(Ze,inf)\n"
		<< L << "    if there are bounds on the problem then\n"
		<< L << "        Find the largest (in magnitude) positive (u_pos) and\n"
		<< L << "        negative (u_neg) step u where the slightly relaxed variable\n"
		<< L << "        bounds:\n"
		<< L << "            xl - delta <= x_k + u * Ze <= xu + delta\n"
		<< L << "        are strictly satisfied (where delta = max_var_bounds_viol).\n"
		<< L << "        if u_pos < u then\n"
		<< L << "            u = u_pos\n"
		<< L << "        end\n"
		<< L << "        if abs(u_neg) < u then\n"
		<< L << "            u = u_neg\n"
		<< L << "        end\n"
		<< L << "    end\n"
		<< L << "    x_fd = x_k + u * Ze\n"
		<< L << "    rGf_fd = ( Z_k' * Gf(x_fd) - rGf_k ) / u\n"
		<< L << "    if norm(rGf_fd,inf) <= min_diag then\n"
		<< L << "        rHL_k = min_diag * eye(n-r)\n"
		<< L << "    else\n"
		<< L << "        if initialization_method == SCALE_IDENTITY then\n"
		<< L << "            rHL_k = norm(rGf_fd,inf) * eye(n-r)\n"
		<< L << "        else if initialization_method == SCALE_DIAGONAL or SCALE_DIAGONAL_ABS then\n"
		<< L << "            *** Make sure that all of the diagonal elements are\n"
		<< L << "            *** positive and that the smallest element is\n"
		<< L << "            *** no smaller than norm(rGf_fd,inf) / max_cond\n"
		<< L << "            *** So that rHL_k will be positive definite an\n"
		<< L << "            *** well conditioned\n"
		<< L << "            min_ele = max( norm(rGf_fd,inf)/max_cond, min_diag )\n"
		<< L << "            if initialization_method == SCALE_DIAGONAL then\n"
		<< L << "                for i = 1 ... n-r\n"
		<< L << "                   diag(i) = max( rGf_fd(i), min_ele )\n"
		<< L << "                end\n"
		<< L << "            else *** SCALE_DIAGONAL_ABS\n"
		<< L << "                for i = 1 ... n-r\n"
		<< L << "                   diag(i) = max( abs(rGf_fd(i)), min_ele )\n"
		<< L << "                end\n"
		<< L << "            end\n"
		<< L << "            rHL_k = diag(diag)\n"
		<< L << "        end\n"
		<< L << "    end\n"
		<< L << "end\n";
}
