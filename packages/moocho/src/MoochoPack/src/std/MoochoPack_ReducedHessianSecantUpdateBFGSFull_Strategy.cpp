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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSFull_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateBFGSFull_Strategy::ReducedHessianSecantUpdateBFGSFull_Strategy(
	const bfgs_update_ptr_t&      bfgs_update
	)
	:bfgs_update_(bfgs_update)
{}

bool ReducedHessianSecantUpdateBFGSFull_Strategy::perform_update(
	VectorWithOpMutable     *s_bfgs
	,VectorWithOpMutable    *y_bfgs
	,bool                   first_update
	,std::ostream           & out
	,EJournalOutputLevel    olevel
	,rSQPAlgo               *algo
	,rSQPState              *s
	,MatrixSymWithOp        *rHL_k
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

}  // end namespace ReducedSpaceSQPPack
