// ////////////////////////////////////////////////////////////////////////////
// PostEvalNewPointBarrier_Step.cpp
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
#include <typeinfo>
#include <iostream>

#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/std/PostEvalNewPointBarrier_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"

#include "StringToBool.h"

#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

bool PostEvalNewPointBarrier_Step::do_step(
  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using DynamicCastHelperPack::dyn_cast;
	using GeneralIterationPack::print_algorithm_step;
	using AbstractLinAlgPack::inv_of_difference;
	using AbstractLinAlgPack::correct_upper_bound_multipliers;
	using AbstractLinAlgPack::correct_lower_bound_multipliers;
	using LinAlgOpPack::Vp_StV;

	rSQPAlgo            &algo   = dyn_cast<rSQPAlgo>(_algo);
	ipState             &s      = dyn_cast<ipState>(_algo.state());
	NLP                 &nlp    = algo.nlp();
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	if(!nlp.is_initialized())
		nlp.initialize(algo.algo_cntr().check_results());

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	IterQuantityAccess<VectorWithOpMutable> &x_iq = s.x();
	IterQuantityAccess<MatrixSymDiagonalStd> &Vl_iq = s.Vl();
	IterQuantityAccess<MatrixSymDiagonalStd> &Vu_iq = s.Vu();

	///***********************************************************
	// Calculate invXl = diag(1/(x-xl)) 
	//  and invXu = diag(1/(xu-x)) matrices
	///***********************************************************

	// get references to x, invXl, and invXu
	VectorWithOpMutable& x = x_iq.get_k(0);
	MatrixSymDiagonalStd& invXu = s.invXu().set_k(0);
	MatrixSymDiagonalStd& invXl = s.invXl().set_k(0);
	
	//std::cout << "xu=\n";
	//nlp.xu().output(std::cout);

	inv_of_difference(1.0, nlp.xu(), x, &invXu.diag());
	inv_of_difference(1.0, x, nlp.xl(), &invXl.diag());

	//std::cout << "invXu=\v";
	//invXu.output(std::cout);

	//std::cout << "\ninvXl=\v";
	//invXl.output(std::cout);
	
	// Check for divide by zeros - 
    using AbstractLinAlgPack::assert_print_nan_inf;
    assert_print_nan_inf(invXu.diag(), "invXu", true, &out); 
    assert_print_nan_inf(invXl.diag(), "invXl", true, &out); 
	// These should never go negative either - could be a useful check

	// Initialize Vu and Vl if necessary
	if ( /*!Vu_iq.updated_k(0) */ Vu_iq.last_updated() == IterQuantity::NONE_UPDATED )
		{
		VectorWithOpMutable& vu = Vu_iq.set_k(0).diag();		
		vu = 0;
		Vp_StV(&vu, s.barrier_parameter().get_k(-1), invXu.diag());

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
			{
			out << "\nInitialize Vu with barrier_parameter * invXu ...\n";
			}
		}

if ( /*!Vl_iq.updated_k(0) */ Vl_iq.last_updated() == IterQuantity::NONE_UPDATED  )
		{
		VectorWithOpMutable& vl = Vl_iq.set_k(0).diag();
		vl = 0;
		Vp_StV(&vl, s.barrier_parameter().get_k(-1), invXl.diag());

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
			{
			out << "\nInitialize Vl with barrier_parameter * invXl ...\n";
			}
		}

	if (s.num_basis().updated_k(0))
		{
		// Basis changed, reorder Vl and Vu
		if (Vu_iq.updated_k(-1))
			{ Vu_iq.set_k(0,-1); }
		if (Vl_iq.updated_k(-1))
			{ Vl_iq.set_k(0,-1); }
			
		VectorWithOpMutable& vu = Vu_iq.set_k(0).diag();
		VectorWithOpMutable& vl = Vl_iq.set_k(0).diag();

		s.P_var_last().permute( BLAS_Cpp::trans, &vu ); // Permute back to original order
		s.P_var_last().permute( BLAS_Cpp::trans, &vl ); // Permute back to original order

		if( (int)olevel >= (int)PRINT_VECTORS ) 
			{
			out	<< "\nx resorted vl and vu to the original order\n" << x;
			}

		s.P_var_current().permute( BLAS_Cpp::no_trans, &vu ); // Permute to new (current) order
		s.P_var_current().permute( BLAS_Cpp::no_trans, &vl ); // Permute to new (current) order

		if( (int)olevel >= (int)PRINT_VECTORS ) 
			{
			out	<< "\nx resorted to new basis \n" << x;
			}
		}

	correct_upper_bound_multipliers(nlp.xu(), +NLP::infinite_bound(), &Vu_iq.get_k(0).diag());
	correct_lower_bound_multipliers(nlp.xl(), -NLP::infinite_bound(), &Vl_iq.get_k(0).diag());
	
	if( (int)olevel >= (int)PRINT_VECTORS ) 
		{
		out << "x=\n"  << s.x().get_k(0);
		out << "xl=\n" << nlp.xl();
		out << "vl=\n" << s.Vl().get_k(0).diag();
		out << "xu=\n" << nlp.xu();
		out << "vu=\n" << s.Vu().get_k(0).diag();
		}

	// Update general algorithm bound multipliers
	VectorWithOpMutable& v = s.nu().set_k(0);
	v = Vu_iq.get_k(0).diag();
	Vp_StV(&v,-1.0,Vl_iq.get_k(0).diag());

	// Print vector information
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out	<< "invXu_k.diag()=\n" << invXu.diag();
		out	<< "invXl_k.diag()=\n" << invXl.diag();
		out	<< "Vu_k.diag()=\n"    << Vu_iq.get_k(0).diag();
		out	<< "Vl_k.diag()=\n"    << Vl_iq.get_k(0).diag();
		out << "nu_k=\n"           << s.nu().get_k(0);
		}

	return true;
	}


void PostEvalNewPointBarrier_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	//const rSQPAlgo   &algo = rsqp_algo(_algo);
	//const rSQPState  &s    = algo.rsqp_state();
	out << L << "# Evaluate information specific to primal / dual barrier algorithms (Post EvalNewPoint)\n"
		<< L << "invXl_k = diag(i, 1/(x(i)-xl))"
		<< L << "invXu_k = diag(i, 1/(xu-x(i)))\n"
		<< L << "if (Vu_k not updated) then\n"
		<< L << "   Vu_k = mu_k*invXu_k\n"
		<< L << "end\n"
		<< L << "if (Vl_k not updated) then\n"
		<< L << "   Vl_k = mu_k*invXl_k\n"
		<< L << "end\n"
		<< L << "nu_k_k = Vu_k.diag() - Vl_k.diag()\n";
	}

} // end namespace ReducedSpaceSQPPack
