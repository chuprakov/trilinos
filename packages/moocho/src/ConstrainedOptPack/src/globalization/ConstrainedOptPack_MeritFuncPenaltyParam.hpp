// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncPenaltyParam.h

#ifndef MERIT_FUNC_PENALTY_PARAM_H
#define MERIT_FUNC_PENALTY_PARAM_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** This class provides interface for setting and retrieving a penalty parameter
  * that many merit functions use {abstract}.
  */
class MeritFuncPenaltyParam {
public:

	///
	virtual ~MeritFuncPenaltyParam() {}

	/// Set the penalty parameter mu
	virtual void mu(value_type mu) = 0;

	/// Get the value of mu
	virtual value_type mu() const = 0;

};	// end class MeritFuncPenaltyParam

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_PENALTY_PARAM_H