// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncPenaltyParams.h

#ifndef MERIT_FUNC_PENALTY_PARAMS_H
#define MERIT_FUNC_PENALTY_PARAMS_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** This class provides interface for setting and retrieving a penalty parameter
  * that many merit functions use {abstract}.
  */
class MeritFuncPenaltyParams {
public:

	///
	class CanNotResize : public std::logic_error
	{public: CanNotResize(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	virtual ~MeritFuncPenaltyParams() {}

	///
	/** Resize the vector of penalty parameters.
	  *
	  * If the subclass can not do this then a #CanNotResize# exception
	  * will be thrown.
	  */
	virtual void resize( size_type n ) = 0;

	/// Get the vector of penalty parameters for setting them
	virtual VectorSlice mu() = 0;

	/// Get the vector of penalty parameters for viewing them
	virtual const VectorSlice mu() const = 0;

};	// end class MeritFuncPenaltyParams

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_PENALTY_PARAMS_H