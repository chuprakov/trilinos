// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_Strategy.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_PROJECTED_FULL_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_PROJECTED_FULL_STRATEGY_H

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/quasi_newton_stats.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Perform BFGS updates on only the free independent (super basic) variables.
 *
 * This method should be very efficient for few super basic variables.
 */
class ReducedHessianSecantUpdateBFGSProjected_Strategy : public ReducedHessianSecantUpdate_Strategy
{
public:
	
	///
	/** <<std comp>> members for the strategy object that will
	 * perform the guts of the BFGS update.
	 */
	STANDARD_COMPOSITION_MEMBERS( BFGSUpdate_Strategy, bfgs_update )

    ReducedHessianSecantUpdateBFGSProjected_Strategy(
		const bfgs_update_ptr_t&      bfgs_update = NULL
		);      

	///
	bool perform_update(
		VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
		,std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,MatrixWithOp *rHL_k
		);
	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

private:
	quasi_newton_stats_iq_member	quasi_newton_stats_;

}; // end class ReducedHessianSecantUpdateBFGSProjected_Strategy

}  // end namespace ReducedSpaceSQPPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_PROJECTED_FULL_STRATEGY_H
