// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_CrossTermExactStd_Step.hpp
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

#ifndef CROSS_TERM_EXACT_STD_STEP_H
#define CROSS_TERM_EXACT_STD_STEP_H

#include "../rSQPAlgo_Step.h"

namespace MoochoPack {

///
/** w_k = Z_k' * HL_k * Ypy_k
  */
class CrossTermExactStd_Step : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CrossTermExactStd_Step

}	// end namespace MoochoPack 

#endif	// CROSS_TERM_EXACT_STD_STEP_H
