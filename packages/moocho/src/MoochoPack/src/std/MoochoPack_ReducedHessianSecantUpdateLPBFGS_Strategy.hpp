// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateLPBFGS_Strategy.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_H

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
class ReducedHessianSecantUpdateLPBFGS_Strategy : public ReducedHessianSecantUpdate_Strategy
{
public:
	
	///
	/** <<std comp>> members for the strategy object that will
	 * perform dense projected BFGS updating.
	 */
	STANDARD_COMPOSITION_MEMBERS( ReducedHessianSecantUpdate_Strategy, dense_proj_update )

	///
	/** <<std comp>> members for the strategy object that will
	 * perform the guts of the BFGS update.
	 */
	STANDARD_COMPOSITION_MEMBERS( BFGSUpdate_Strategy, bfgs_update )

	///
	/** Set the ratio of the number of active independent varaibles in the
	  * active-set of the last two calls before limited memory projected updating
	  * for superbasic variables only is started.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, act_set_frac_proj_start )

	///
	/** Set the tolerance for Langrange multipliers for fixed variables
	 * below which rows/cols from rHL_RR will not be dropped.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, super_basic_mult_drop_tol )
		
	///
	/** Set the maximum number of super basic variables under which switching
	 * to dense projected PBFGS updating will be allowed. 
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_superbasics_switch_dense )

	///
	/** Set the ratio of the number of inequality constraints in the
	  * active-set of the last two calls before considering switching to
	  * dense projected updating for superbasic variables.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, act_set_frac_switch_dense )

	///
	/** Set the minimum number of BFGS updates to perform on the LBFGS matrix
	 * before considering switching to the dense projected BFGS updating.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, min_num_updates_switch_dense )

	///
	/** Set the maximum number of BFGS updates to perform on the LBFGS matrix
	 * before automatically switching to the dense projected BFGS updating
	 * reguardless if the active set has calmed down or not but only if
	 * the number of fixed independent variables is sufficiently large.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, max_num_updates_switch_dense )

	///
	/** Set maximum number of previous BFGS updates to initialize the new dense
	 * projected BFGS matrix with.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_add_recent_updates )

	///
    ReducedHessianSecantUpdateLPBFGS_Strategy(
		const dense_proj_update_ptr_t&  dense_proj_update            = NULL
		,const bfgs_update_ptr_t&       bfgs_update                  = NULL
		,value_type                     act_set_frac_proj_start      = 0.5
		,value_type                     super_basic_mult_drop_tol    = 1e-5
		,size_type                      num_superbasics_switch_dense = 500
		,value_type                     act_set_frac_switch_dense    = 0.8
		,size_type                      min_num_updates_switch_dense = 0
		,size_type                      max_num_updates_switch_dense = 999999
		,size_type                      num_add_recent_updates       = 10
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

	// /////////////////////////////
	// Private data members

	quasi_newton_stats_iq_member    quasi_newton_stats_;
	act_set_stats_iq_member         act_set_stats_;

}; // end class ReducedHessianSecantUpdateLPBFGS_Strategy

}  // end namespace ReducedSpaceSQPPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_H
