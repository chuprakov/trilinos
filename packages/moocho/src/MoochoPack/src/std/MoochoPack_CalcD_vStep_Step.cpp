// ////////////////////////////////////////////////////////////////////////////
// CalcD_vStep_Step.cpp
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

#include <limits>
#include <ostream>
#include <iostream>

#include "ReducedSpaceSQPPack/include/std/CalcD_vStep_Step.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
//#include "ConstrainedOptimizationPack/include/print_vector_change_stats.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "dynamic_cast_verbose.h"


bool ReducedSpaceSQPPack::CalcD_vStep_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
	{
	using DynamicCastHelperPack::dyn_cast;
	using GeneralIterationPack::print_algorithm_step;
	using AbstractLinAlgPack::ele_wise_prod;

	rSQPAlgo &algo = rsqp_algo(_algo);
	ipState	&s = dyn_cast<ipState>(_algo.state());
	NLPFirstOrderInfo   &nlp    = dyn_cast<NLPFirstOrderInfo>(algo.nlp());

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
		}

	// Get iteration quantities
	const value_type& mu = s.barrier_parameter().get_k(0);
	const VectorWithOp &d_k = s.d().get_k(0);
	const MatrixSymDiagonalStd& invXl = s.invXl().get_k(0);
	const MatrixSymDiagonalStd& invXu = s.invXu().get_k(0);
	const MatrixSymDiagonalStd& Vl = s.Vl().get_k(0);
	const MatrixSymDiagonalStd& Vu = s.Vu().get_k(0);

	VectorWithOpMutable& dvl_k = s.dvl().set_k(0);
	VectorWithOpMutable& dvu_k = s.dvu().set_k(0);

	lowerbound_multipliers_step(mu, invXl.diag(), Vl.diag(), d_k, &dvl_k);
	upperbound_multipliers_step(mu, invXu.diag(), Vu.diag(), d_k, &dvu_k);

	/*
	// dvl = mu*invXl*e - vl - invXl*Vl*d_k
	dvl_k = 0;
	ele_wise_prod(-1.0, invXl.diag(), Vl.diag(), &dvl_k);
	ele_wise_prod(1.0, dvl_k, d_k, &dvl_k);

	std::cout << "d_k =\n" << d_k;
 	std::cout << "-invXl*Vl*d_k = \n" << dvl_k;
 
	dvl_k.axpy(-1.0, Vl.diag());
	
 	std::cout << "-vl-invXl*Vl*d_k = \n" << dvl_k;

	dvl_k.axpy(mu, invXl.diag());

	 std::cout << "dvl_k = \n" << dvl_k;

	// dvu = mu*invXu*e - vu + invXu*Vu*d_k
	dvu_k = 0;
	ele_wise_prod(1.0, invXu.diag(), Vu.diag(), &dvu_k);
	ele_wise_prod(1.0, dvu_k, d_k, &dvu_k);

	dvu_k.axpy(-1.0, Vu.diag());
	
	dvu_k.axpy(mu, invXu.diag());
	*/
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out	<< "\nx_k = \n" << s.x().get_k(0)
			<< "\nxl = \n" << nlp.xl()
			<< "\nxu = \n" << nlp.xu()
			<< "\ndvl_k = \n" << dvl_k
			<< "\ndvu_k = \n" << dvu_k;
		}

	return true;
	}

void ReducedSpaceSQPPack::CalcD_vStep_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculates the search direction for the dual variables\n"
		<< L << "dvl_k = mu*invXl_k*e - vl_k - invXl_k*Vl_k*d_k\n"
		<< L << "dvu_k = mu*invXu_k*e - vu_k + invXu_k*Vu_k*d_k\n";
}
