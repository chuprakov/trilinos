// ////////////////////////////////////////////////////////////////////////////
// LineSearchWatchDog_Step.h

#ifndef LINE_SEARCH_WATCH_DOG_STEP_H
#define LINE_SEARCH_WATCH_DOG_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearch_Strategy.h"
#include "LinAlgPack/include/VectorClass.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardAggregationMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Implements watchdog line search.
  */
class LineSearchWatchDog_Step : public LineSearch_Step {
public:

	/// <<std comp>> members for direct_line_search
	STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search)

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	LineSearchWatchDog_Step(
			  const direct_line_search_ptr_t&	direct_line_search	= 0
			, const merit_func_ptr_t&			merit_func			= 0
			, value_type						use_line_search_correct_kkt_tol = 0.0
			, value_type						eta					= 1.0e-4	);

	/// Set the KKT tolerance below which we will start applying the watchdog
	void set_use_line_search_correct_kkt_tol( value_type use_line_search_correct_kkt_tol )
	{	use_line_search_correct_kkt_tol_ = use_line_search_correct_kkt_tol; }
	///
	value_type use_line_search_correct_kkt_tol() const
	{	return use_line_search_correct_kkt_tol_; }

	/// Set the Armijo test parameter eta
	void set_eta(value_type eta)
	{	eta_ = eta; }
	///
	value_type eta() const
	{	return eta_; }

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	value_type			use_line_search_correct_kkt_tol_;
	value_type			eta_;
	int					watch_k_;	// >= 0 means that we are using the watchdog.
	Vector				xo_;
	value_type			fo_;
	value_type			nrm_co_;
	Vector				do_;
	value_type			phio_;
	value_type			Dphio_;
	value_type			phiop1_;
	
};	// end class LineSearchWatchDog_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// LINE_SEARCH_WATCH_DOG_STEP_H