// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalcNLP.h

#ifndef MERIT_FUNC_CALC_NLP_H
#define MERIT_FUNC_CALC_NLP_H

#include "MeritFuncCalc.h"
#include "MeritFuncNLP.h"
#include "Misc/include/StandardAggregationMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Adds the ability to compute phi(f(x), c(x)) at x
  * directly instead of having to compute f, and c first.
  * This class uses an aggregate NLP to perform the computations of f(x)
  * and c(x).
  */
class MeritFuncCalcNLP : public MeritFuncCalc {
public:

	/// <<std aggr>> stereotype members for phi.
	STANDARD_CONST_AGGREGATION_MEMBERS( MeritFuncNLP, phi )

	/// <<std aggr>> stereotype members for nlp.
	STANDARD_CONST_AGGREGATION_MEMBERS( NLP, nlp )

	///
	MeritFuncCalcNLP( const MeritFuncNLP* phi = 0, const NLP* nlp = 0 );

	// ////////////////////////////////////////////
	// Overridden from MeritFuncCalc

	///
	/** Return the value of the merit function at x.
	  * Here phi(x) is calculated directly using the nlp.
	  */
	value_type operator()(const VectorSlice& x) const;

	/// Calls phi().deriv() on phi.
	value_type deriv() const;

	/// Calls phi().print_merit_func(....).
	void print_merit_func(std::ostream& out
		, const std::string& leading_str) const;

};	// end class MeritFuncCalcNLP

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_CALC_NLP_H