// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_NumFixedDepIndep_AddedStep.hpp
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

#ifndef NUM_FIXED_DEP_INDEP_ADDED_STEP_H
#define NUM_FIXED_DEP_INDEP_ADDED_STEP_H

#include "../rSQPAlgo_Step.h"

namespace MoochoPack {

///
/** Computes and outputs the number of fixed variables from the dependent
  * and independent set..
  *
  * These statistics only make sense for variable reduction decompositions.
  */
class NumFixedDepIndep_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class NumFixedDepIndep_AddedStep

}	// end namespace MoochoPack 

#endif	// NUM_FIXED_DEP_INDEP_ADDED_STEP_H
