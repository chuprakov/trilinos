// ////////////////////////////////////////////////////////////////////////////
// LineSearchFullStepAfterKIter_Step.h

#ifndef LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
#define LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H

#include <limits>

#include "../rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearch_Strategy.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardAggregationMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Changes from a line search step to just taking full steps after
  * full_steps_after_k iterations.
  */
class LineSearchFullStepAfterKIter_Step : public LineSearch_Step {
public:

	/// <<std comp>> members for the line search step
	STANDARD_COMPOSITION_MEMBERS(LineSearch_Step,line_search)

	///
	LineSearchFullStepAfterKIter_Step(
			  const line_search_ptr_t&	line_search			= 0
			, int						full_steps_after_k
											= std::numeric_limits<int>::max()	)
		: line_search_(line_search)
			, full_steps_after_k_(full_steps_after_k)
	{}

	/// 
	void full_steps_after_k( int full_steps_after_k )
	{	full_steps_after_k_ = full_steps_after_k; }
	///
	value_type full_steps_after_k() const
	{	return full_steps_after_k_; }

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	int		full_steps_after_k_;
	
};	// end class LineSearchFullStepAfterKIter_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H