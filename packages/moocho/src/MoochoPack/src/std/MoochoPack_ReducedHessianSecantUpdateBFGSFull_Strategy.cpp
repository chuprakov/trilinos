// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSFull_Strategy.cpp
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

#include "MoochoPack_ReducedHessianSecantUpdateBFGSFull_Strategy.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

ReducedHessianSecantUpdateBFGSFull_Strategy::ReducedHessianSecantUpdateBFGSFull_Strategy(
	const bfgs_update_ptr_t&      bfgs_update
	)
	:bfgs_update_(bfgs_update)
{}

bool ReducedHessianSecantUpdateBFGSFull_Strategy::perform_update(
	VectorMutable           *s_bfgs
	,VectorMutable          *y_bfgs
	,bool                   first_update
	,std::ostream           & out
	,EJournalOutputLevel    olevel
	,NLPAlgo                *algo
	,NLPAlgoState           *s
	,MatrixSymOp            *rHL_k
	)
{
	bfgs_update().perform_update(
		s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
		,rHL_k, &quasi_newton_stats_(*s).set_k(0)
		);
	return true;
}

void ReducedHessianSecantUpdateBFGSFull_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Perform BFGS update on full matrix where: B = rHL_k\n";
	bfgs_update().print_step(out,L);
}

}  // end namespace MoochoPack
