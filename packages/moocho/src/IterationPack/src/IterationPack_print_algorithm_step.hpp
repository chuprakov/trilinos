// //////////////////////////////////////////////////////////////////
// print_algorithm_step.h

#ifndef PRINT_ALGORITHM_STEP_H
#define PRINT_ALGORITHM_STEP_H

#include <iosfwd>

#include "Algorithm.h"

namespace GeneralIterationPack {

///
/** Prints to 'out' the algorithm step.
  *
  * 
  *
  */
void print_algorithm_step( const Algorithm& algo, Algorithm::poss_type step_poss
	, EDoStepType type, Algorithm::poss_type assoc_step_poss
	, std::ostream& out );

}	// end namespace GeneralIterationPack

#endif