// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceIP_Strategy.cpp
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
//#include <limits>
//#include <sstream>

#include "ReducedSpaceSQPPack/include/std/CheckConvergenceIP_Strategy.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

CheckConvergenceIP_Strategy::CheckConvergenceIP_Strategy(
	EOptErrorCheck         opt_error_check
	,EScaleKKTErrorBy      scale_opt_error_by
	,EScaleKKTErrorBy      scale_feas_error_by
	,EScaleKKTErrorBy      scale_comp_error_by
	,bool                  scale_opt_error_by_Gf
	)
	:
	CheckConvergenceStd_Strategy(
	  opt_error_check,
	  scale_opt_error_by,
	  scale_feas_error_by,
	  scale_comp_error_by,
	  scale_opt_error_by_Gf
	  )
	{}

bool CheckConvergenceIP_Strategy::Converged(
  Algorithm& _algo
  )
	{
	using DynamicCastHelperPack::dyn_cast;

	// Calculate kkt errors and check for overall convergence
	//bool found_solution = CheckConvergenceStd_Strategy::Converged(_algo);
	bool found_solution = false;

	// Recalculate the complementarity error including mu
	
	// Get the iteration quantities
	ipState &s = dyn_cast<ipState>(*_algo.get_state());
	rSQPAlgo& algo = rsqp_algo(_algo);
	NLP& nlp = algo.nlp();
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// Get necessary iteration quantities
	const value_type &mu_km1 = s.barrier_parameter().get_k(-1);
	const VectorWithOp& x_k = s.x().get_k(0);
	const VectorWithOpMutable& Gf_k = s.Gf().get_k(0);
	const VectorWithOp& rGL_k = s.rGL().get_k(0);
	const VectorWithOp& c_k = s.c().get_k(0);
	const VectorWithOp& vl_k = s.Vl().get_k(0).diag();
	const VectorWithOp& vu_k = s.Vu().get_k(0).diag();
	
	// Calculate the errors with Andreas' scaling
	value_type& opt_err = s.opt_kkt_err().set_k(0);
	value_type& feas_err = s.feas_kkt_err().set_k(0);
	value_type& comp_err = s.comp_kkt_err().set_k(0);

	// scaling
	value_type scale_1 = 1 + x_k.norm_1()/x_k.dim();

	MemMngPack::ref_count_ptr<VectorWithOpMutable> temp = Gf_k.clone();
	temp->axpy(-1.0, vl_k);
	temp->axpy(1.0, vu_k);
	value_type scale_2 = temp->norm_1();
	scale_2 += vl_k.norm_1() + vu_k.norm_1();

	*temp = nlp.infinite_bound();
	const size_type nlb = num_bounded(nlp.xl(), *temp, nlp.infinite_bound());
	*temp = -nlp.infinite_bound();
	const size_type nub = num_bounded(*temp, nlp.xu(), nlp.infinite_bound());
	scale_2 = 1 + scale_2/(1+nlp.m()+nlb+nub);

	// Calculate the opt_err
	opt_err = rGL_k.norm_inf() / scale_2;

	// Calculate the feas_err
	feas_err = c_k.norm_inf() / scale_1;
	
	// Calculate the comp_err
	comp_err = 0.0;

	out << "x=\n";
	s.x().get_k(0).output(out);
	out << "xl=\n";
	nlp.xl().output(out);
	out << "vl=\n";
	s.Vl().get_k(0).diag().output(out);
	out << "xu=\n";
	nlp.xu().output(out);
	out << "vu=\n";
	s.Vu().get_k(0).diag().output(out);

	comp_err = IP_comp_err_with_mu(mu_km1, nlp.infinite_bound(), s.x().get_k(0), nlp.xl(), nlp.xu(), s.Vl().get_k(0).diag(), s.Vu().get_k(0).diag());
	comp_err = comp_err / scale_2;

	// check for convergence
	
	const value_type opt_tol = algo.algo_cntr().opt_tol();
	const value_type feas_tol = algo.algo_cntr().feas_tol();
	const value_type comp_tol = algo.algo_cntr().comp_tol();
	if (opt_err < opt_tol && feas_err < feas_tol && comp_err < comp_tol && mu_km1 < comp_tol)
		{
		found_solution = true;
		}

	if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) || (int(olevel) >= int(PRINT_BASIC_ALGORITHM_INFO) && found_solution) )
		{
		out	
			<< "\nopt_kkt_err_k   = " << opt_err << ( opt_err < opt_tol ? " < " : " > " )
			<< "opt_tol = " << opt_tol
			<< "\nfeas_kkt_err_k   = " << feas_err << ( feas_err < feas_tol ? " < " : " > " )
			<< "feas_tol = " << feas_tol
			<< "\ncomp_kkt_err_k   = " << comp_err << ( comp_err < comp_tol ? " < " : " > " )
			<< "comp_tol = " << comp_tol
			<< "\nbarrier_parameter_k (mu_km1) = " << mu_km1 << ( mu_km1 < comp_tol ? " < " : " > " )
			<< "comp_tol = " << comp_tol
			<< std::endl;
		}
		
	return found_solution;
	}

void CheckConvergenceIP_Strategy::print_step( const Algorithm& _algo, std::ostream& out, const std::string& L ) const
	{
	out 
		<< L << "CheckConvergenceIP_Strategy\n"
		<< L << " IP_found_solution = CheckConvergedStd_Strategy::Converged(_algo, reportFinalSolution)\n";
	
	CheckConvergenceStd_Strategy::print_step(_algo, out, L+"   ");
	
	out 
		<< L << "# recalculate comp_err\n"
		<< L << "comp_err_k = 0.0"
		<< L << "for all i = 1 to n\n"
		<< L << "   comp_err_k = max( comp_err_k, vl_k(i)*(x_k(i)-xl_k(i))-mu_km1, vu_k(i)*(xu_k(i)-x_k(i))-mu_k )\n"
		<< L << "next i\n"
		<< L << "if IP_found_solution then\n"
		<< L << "   IP_found_solution = false\n"
		<< L << "   if (comp_err_k < comp_tol && mu_k < comp_tol) then\n"
		<< L << "      IP_found_solution = true\n"
		<< L << "   end\n"
		<< L << "end\n"
		<< L << "return IP_found_solution\n";
	}

}	// end namespace ReducedSpaceSQPPack

