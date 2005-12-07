// //////////////////////////////////////////////////////////////////
// IterationPack_print_algorithm_step.hpp
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

#ifndef PRINT_ALGORITHM_STEP_H
#define PRINT_ALGORITHM_STEP_H

#include <iosfwd>

#include "IterationPack_Algorithm.hpp"

namespace IterationPack {

///
/** Prints to 'out' the algorithm step.
  *
  * This function can be used by \c AlgorithmStep subclasses in their
  * \c AlgorithmStep::print_step() methods to print the header for the
  * algorithm step.  This will print the number of the step in the algorithm,
  * its name, and the name of its concreate subclass.
  */
void print_algorithm_step( const Algorithm& algo, Algorithm::poss_type step_poss
	, EDoStepType type, Algorithm::poss_type assoc_step_poss
	, std::ostream& out );

}	// end namespace IterationPack

#endif // PRINT_ALGORITHM_STEP_H
