// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_Strategy.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_PROJECTED_FULL_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_PROJECTED_FULL_STRATEGY_H

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/quasi_newton_stats.h"
#include "ReducedSpaceSQPPack/include/std/act_set_stats.h"
#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasic.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

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

	///
	/** Set the ratio of the number of inequality constraints in the
	  * active-set of the last two calls before a projected updating
	  * for superbasic variables only is started.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, proj_start_frac )

	///
	/** Set the tolerance for Langrange multipliers for fixed variables
	 * below which rows/cols from rHL_RR will not be dropped.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, super_basic_mult_drop_tol )
		
	///
    ReducedHessianSecantUpdateBFGSProjected_Strategy(
		const bfgs_update_ptr_t&      bfgs_update                = NULL
		,value_type                   proj_start_frac            = 0.8
		,value_type                   super_basic_mult_drop_tol  = 1e-5
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
	
	// //////////////////////////////
	// Private types

	typedef std::vector<size_type>                         i_x_t;
	typedef std::vector<MatrixHessianSuperBasic::EBounds>  bnd_fixed_t;

	// /////////////////////////////
	// Private data members

	i_x_t                           i_x_free_;   // Keeps track of who is free
	i_x_t                           i_x_fixed_;  // Keeps track of who is fixed
	bnd_fixed_t                     bnd_fixed_;
	quasi_newton_stats_iq_member    quasi_newton_stats_;
	act_set_stats_iq_member         act_set_stats_;

}; // end class ReducedHessianSecantUpdateBFGSProjected_Strategy

}  // end namespace ReducedSpaceSQPPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_PROJECTED_FULL_STRATEGY_H
