// ////////////////////////////////////////////////////////
// ComputeMinMult.h

#ifndef COMPUTE_MIN_MULT_H
#define COMPUTE_MIN_MULT_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

/** @name Compute the minimum absolute value of the given
  * Lagrange multipliers.
  *
  * A small Lagrange multiplier indicates degeneracy.
  */
//@{

/// Minimum |mu(i)|
value_type min_abs( const VectorSlice& mu );

/// Minimum |mu(i)|
value_type min_abs( const SpVectorSlice& mu );

//@}

}	// end namespace ConstrainedOptimizationPack

#endif	// COMPUTE_MIN_MULT_H
