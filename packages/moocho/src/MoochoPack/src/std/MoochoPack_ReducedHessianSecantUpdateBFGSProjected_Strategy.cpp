// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_Strategy.cpp

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSProjected_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateBFGSProjected_Strategy::ReducedHessianSecantUpdateBFGSProjected_Strategy(
	const bfgs_update_ptr_t&      bfgs_update
	)
	:
	    bfgs_update_(bfgs_update)
{}

bool ReducedHessianSecantUpdateBFGSProjected_Strategy::perform_update(
	VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
	,std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,MatrixWithOp *rHL_k
	)
{
	assert(0); // Implement this!

	return true;
}

void ReducedHessianSecantUpdateBFGSProjected_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Perform BFGS update on only free independent (super basic) variables.\n"
		<< L << "ToDo: Finish implementation!\n";
}

}  // end namespace ReducedSpaceSQPPack
