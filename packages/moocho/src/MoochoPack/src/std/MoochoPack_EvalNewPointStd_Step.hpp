// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointStd_Step.h

#ifndef EVAL_NEW_POINT_STD_STEP_H
#define EVAL_NEW_POINT_STD_STEP_H

#include "../rSQPAlgo_StepBaseClasses.h"

namespace ReducedSpaceSQPPack {

///
/** Standard new point evaluation step class.
  *
  * This class calcualtes Gc, updates Z, and Y, and calculates Gf, c, and f in that order.
  * Also the lagrange multipliers (lambda) for the equality constraints (c) can be
  * calculated if compute_lambda(true) is called.
  */
class EvalNewPointStd_Step : public EvalNewPoint_Step {
public:

	/// set new_point == true by default.
	EvalNewPointStd_Step()
		: new_point_(true)
	{}

	///
	/** Call to set #newx# = #new_point# which is passed to the to the first nlp calc with is Gc.
	  *
	  * This is primarily used when a new basis must be selected for variable
	  * reduction decompositions for Z, and Y.
	  *
	  * For this operation to have its effect it must be called before each call to
	  * #do_step()# and therefore the effect does not presist between calls of
	  * #do_step()#.
	  */
	void set_new_point(bool new_point) {
		new_point_ = new_point;
	}

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
	bool		new_point_;	// flag for if a newx will be calculated for
	
};	// end class EvalNewPointStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// EVAL_NEW_POINT_STD_STEP_H