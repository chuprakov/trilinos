// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_Strategy.cpp
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

#include <assert.h>

#include <ostream>
#include <limits>
#include <sstream>

#include "ReducedSpaceSQPPack/src/std/CheckConvergenceStd_Strategy.hpp"
#include "ReducedSpaceSQPPack/src/rSQPAlgoContainer.hpp"
#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackExceptions.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"
#include "AbstractLinAlgPack/src/MatrixOpOut.hpp"
#include "AbstractLinAlgPack/src/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/VectorOut.hpp"
#include "AbstractLinAlgPack/src/assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ThrowException.hpp"

namespace ReducedSpaceSQPPack {

CheckConvergenceStd_Strategy::CheckConvergenceStd_Strategy(
	EOptErrorCheck         opt_error_check
	,EScaleKKTErrorBy      scale_opt_error_by
	,EScaleKKTErrorBy      scale_feas_error_by
	,EScaleKKTErrorBy      scale_comp_error_by
	,bool                  scale_opt_error_by_Gf
	)
	:
	CheckConvergence_Strategy(
	  opt_error_check,
	  scale_opt_error_by,
	  scale_feas_error_by,
	  scale_comp_error_by,
	  scale_opt_error_by_Gf
	  )
	{}

bool CheckConvergenceStd_Strategy::Converged(
  Algorithm& _algo
  )
	{
	using AbstractLinAlgPack::assert_print_nan_inf;
	using AbstractLinAlgPack::combined_nu_comp_err;
	
	rSQPAlgo	&algo	  = rsqp_algo(_algo);
	rSQPState	&s		  = algo.rsqp_state();
	NLP			&nlp	  = algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	const size_type
		n  = nlp.n(),
		m  = nlp.m(),
		mI = nlp.mI(),
		nb = nlp.num_bounded_x();

	// Get the iteration quantities
	IterQuantityAccess<value_type>
		&opt_kkt_err_iq  = s.opt_kkt_err(),
		&feas_kkt_err_iq = s.feas_kkt_err(),
	    &comp_kkt_err_iq = s.comp_kkt_err();
	
	IterQuantityAccess<VectorMutable>
		&x_iq       = s.x(),
		&d_iq       = s.d(),
		&Gf_iq      = s.Gf(),
		*c_iq       = m  ? &s.c() : NULL,
		*h_iq       = mI ? &s.h() : NULL,
		&rGL_iq     = s.rGL(),
		&GL_iq      = s.GL(),
		*lambda_iq  = m  ? &s.lambda()  : NULL,
		*lambdaI_iq = mI ? &s.lambdaI() : NULL,
		&nu_iq      = s.nu();

	// opt_err = (||rGL||inf or ||GL||) / (||Gf|| + scale_kkt_factor)
	value_type 
		norm_inf_Gf_k = 0.0,
		norm_inf_GLrGL_k = 0.0;

	if( scale_opt_error_by_Gf() ) 
		{
		assert_print_nan_inf( norm_inf_Gf_k = Gf_iq.get_k(0).norm_inf(),
							  "||Gf_k||inf",true,&out);
		}

	// NOTE:
	// The strategy object CheckConvergenceIP_Strategy assumes
	//  that this will always be the gradient of the lagrangian
	//  of the original problem, not the gradient of the lagrangian
	//  for psi. (don't use augmented nlp info here)
	if( opt_error_check() == OPT_ERROR_REDUCED_GRADIENT_LAGR ) 
		{
		assert_print_nan_inf( norm_inf_GLrGL_k = rGL_iq.get_k(0).norm_inf(),
							  "||rGL_k||inf",true,&out);
		}
	else 
		{
		assert_print_nan_inf( norm_inf_GLrGL_k = GL_iq.get_k(0).norm_inf(),
							  "||GL_k||inf",true,&out);
		}

	const value_type
		opt_scale_factor = 1.0 + norm_inf_Gf_k,
		opt_err = norm_inf_GLrGL_k / opt_scale_factor;
	
	// feas_err
	const value_type feas_err = ( ( m ? c_iq->get_k(0).norm_inf() : 0.0 ) );

	// comp_err
	value_type comp_err = 0.0;
	if (nb > 0)
		{
		comp_err = combined_nu_comp_err(nu_iq.get_k(0), x_iq.get_k(0), nlp.xl(), nlp.xu());
		}

	if(m)
		assert_print_nan_inf( feas_err,"||c_k||inf",true,&out);

	// scaling factors
	const value_type 
		scale_opt_factor = CalculateScalingFactor(s, scale_opt_error_by()),
		scale_feas_factor = CalculateScalingFactor(s, scale_feas_error_by()),
		scale_comp_factor = CalculateScalingFactor(s, scale_comp_error_by());

	// kkt_err
	const value_type
		opt_kkt_err_k  = opt_err/scale_opt_factor,
 		feas_kkt_err_k = feas_err/scale_feas_factor,
		comp_kkt_err_k = comp_err/scale_comp_factor;

	// update the iteration quantities
	opt_kkt_err_iq.set_k(0) = opt_kkt_err_k;
	feas_kkt_err_iq.set_k(0) = feas_kkt_err_k;
	comp_kkt_err_iq.set_k(0) = comp_kkt_err_k;

	// step_err
	value_type step_err = 0.0;
	if( d_iq.updated_k(0) ) 
		{
			step_err = AbstractLinAlgPack::max_rel_step(x_iq.get_k(0),d_iq.get_k(0));
			assert_print_nan_inf( step_err,"max(d(i)/max(1,x(i)),i=1...n)",true,&out);
		}

	const value_type
		opt_tol		= algo.algo_cntr().opt_tol(),
		feas_tol	= algo.algo_cntr().feas_tol(),
		comp_tol	= algo.algo_cntr().comp_tol(),
		step_tol	= algo.algo_cntr().step_tol();

	const bool found_solution = 
		opt_kkt_err_k < opt_tol 
		&& feas_kkt_err_k < feas_tol 
		&& comp_kkt_err_k < comp_tol 
		&& step_err < step_tol;

	if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) || (int(olevel) >= int(PRINT_BASIC_ALGORITHM_INFO) && found_solution) )
		{
		out	
			<< "\nscale_opt_factor = " << scale_opt_factor
			<< " (scale_opt_error_by = " << (scale_opt_error_by()==SCALE_BY_ONE ? "SCALE_BY_ONE"
											 : (scale_opt_error_by()==SCALE_BY_NORM_INF_X ? "SCALE_BY_NORM_INF_X"
												: "SCALE_BY_NORM_2_X" ) ) << ")"

			<< "\nscale_feas_factor = " << scale_feas_factor
			<< " (scale_feas_error_by = " << (scale_feas_error_by()==SCALE_BY_ONE ? "SCALE_BY_ONE"
											 : (scale_feas_error_by()==SCALE_BY_NORM_INF_X ? "SCALE_BY_NORM_INF_X"
												: "SCALE_BY_NORM_2_X" ) ) << ")"

			<< "\nscale_comp_factor = " << scale_comp_factor
			<< " (scale_comp_error_by = " << (scale_comp_error_by()==SCALE_BY_ONE ? "SCALE_BY_ONE"
											 : (scale_comp_error_by()==SCALE_BY_NORM_INF_X ? "SCALE_BY_NORM_INF_X"
												: "SCALE_BY_NORM_2_X" ) ) << ")"
			<< "\nopt_scale_factor = " << opt_scale_factor
			<< " (scale_opt_error_by_Gf = " << (scale_opt_error_by_Gf()?"true":"false") << ")"
			<< "\nopt_kkt_err_k    = " << opt_kkt_err_k << ( opt_kkt_err_k < opt_tol ? " < " : " > " )
			<< "opt_tol  = " << opt_tol
			<< "\nfeas_kkt_err_k   = " << feas_kkt_err_k << ( feas_kkt_err_k < feas_tol ? " < " : " > " )
			<< "feas_tol = " << feas_tol
			<< "\ncomp_kkt_err_k   = " << comp_kkt_err_k << ( comp_kkt_err_k < comp_tol ? " < " : " > " )
			<< "comp_tol = " << comp_tol
			<< "\nstep_err         = " << step_err << ( step_err < step_tol ? " < " : " > " )
			<< "step_tol = " << step_tol
			<< std::endl;
		}
	
	return found_solution;

	}

void CheckConvergenceStd_Strategy::print_step( const Algorithm& _algo, std::ostream& out, const std::string& L ) const
	{
	out
		<< L << "*** Check to see if the KKT error is small enough for convergence\n"
		<< L << "if scale_(opt|feas|comp)_error_by == SCALE_BY_ONE then\n"
		<< L << "    scale_(opt|feas|comp)_factor = 1.0\n"
		<< L << "else if scale_(opt|feas|comp)_error_by == SCALE_BY_NORM_2_X then\n"
		<< L << "    scale_(opt|feas|comp)_factor = 1.0 + norm_2(x_k)\n"
		<< L << "else if scale_(opt|feas|comp)_error_by == SCALE_BY_NORM_INF_X then\n"
		<< L << "    scale_(opt|feas|comp)_factor = 1.0 + norm_inf(x_k)\n"
		<< L << "end\n"
		<< L << "if scale_opt_error_by_Gf == true then\n"
		<< L << "    opt_scale_factor = 1.0 + norm_inf(Gf_k)\n"
		<< L << "else\n"
		<< L << "    opt_scale_factor = 1.0\n"
		<< L << "end\n";
	if( opt_error_check() == OPT_ERROR_REDUCED_GRADIENT_LAGR )
		{
		out
			<< L << "opt_err = norm_inf(rGL_k)/opt_scale_factor\n";
		}
	else
		{
		out
			<< L << "opt_err = norm_inf(GL_k)/opt_scale_factor\n";
		}

	out
		<< L << "feas_err = norm_inf(c_k)\n"
		<< L << "comp_err = max(i, nu(i)*(xu(i)-x(i)), -nu(i)*(x(i)-xl(i)))\n"
		<< L << "opt_kkt_err_k = opt_err/scale_opt_factor\n"
		<< L << "feas_kkt_err_k = feas_err/scale_feas_factor\n"
		<< L << "comp_kkt_err_k = feas_err/scale_comp_factor\n"
		<< L << "if d_k is updated then\n"
		<< L << "    step_err = max( |d_k(i)|/(1+|x_k(i)|), i=1..n )\n"
		<< L << "else\n"
		<< L << "    step_err = 0\n"
		<< L << "end\n"
		<< L << "if opt_kkt_err_k < opt_tol\n"
		<< L << "       and feas_kkt_err_k < feas_tol\n"
		<< L << "       and step_err < step_tol then\n"
		<< L << "   report optimal x_k, lambda_k and nu_k to the nlp\n"
		<< L << "   terminate, the solution has beed found!\n"
		<< L << "end\n";
	}


value_type CheckConvergenceStd_Strategy::CalculateScalingFactor( rSQPState& state, EScaleKKTErrorBy scale_by ) const
	{
	// scale_kkt_factor
	value_type scale_factor = 1.0;
	switch(scale_by) 
		{
		case SCALE_BY_ONE:
			scale_factor = 1.0;
			break;
		case SCALE_BY_NORM_2_X:
			scale_factor = 1.0 + state.x().get_k(0).norm_2();
			break;
		case SCALE_BY_NORM_INF_X:
			scale_factor = 1.0 + state.x().get_k(0).norm_inf();
			break;
		default:
			assert(0);	// Should never be called
		}

	return scale_factor;
	}

}	// end namespace ReducedSpaceSQPPack

