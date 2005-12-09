// ////////////////////////////////////////////////////////
// ConstrainedOptPack_ComputeMinMult.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef COMPUTE_MIN_MULT_H
#define COMPUTE_MIN_MULT_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** @name Compute the minimum absolute value of the given
  * Lagrange multipliers.
  *
  * A small Lagrange multiplier indicates degeneracy.
  */
//@{

/// Minimum |mu(i)|
value_type min_abs( const DVectorSlice& mu );

/// Minimum |mu(i)|
value_type min_abs( const SpVectorSlice& mu );

//@}

}	// end namespace ConstrainedOptPack

#endif	// COMPUTE_MIN_MULT_H
