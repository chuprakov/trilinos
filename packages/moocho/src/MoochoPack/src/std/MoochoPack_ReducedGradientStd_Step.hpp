// ////////////////////////////////////////////////////////////////////////////
// ReducedGradientStd_Step.h
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

#ifndef REDUCED_GRADIENT_STD_STEP_H
#define REDUCED_GRADIENT_STD_STEP_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/AlgorithmStep.h"

namespace ReducedSpaceSQPPack {

///
/** Computes the reducecd gradient of the objective <tt>rGf_k = Z_k' * Gf_k</tt>
  */
class ReducedGradientStd_Step
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

};	// end class ReducedGradientStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// REDUCED_GRADIENT_STD_STEP_H
