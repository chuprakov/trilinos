// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSFull_Strategy.cpp

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSFull_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateBFGSFull_Strategy::ReducedHessianSecantUpdateBFGSFull_Strategy(
	const bfgs_update_ptr_t&      bfgs_update
	)
	:
	    bfgs_update_(bfgs_update)
{}

bool ReducedHessianSecantUpdateBFGSFull_Strategy::perform_update(
	VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
	,std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,MatrixWithOp *rHL_k
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
