// ////////////////////////////////////////////////////////////////////////////
// RangeSpaceStepStd_Step.h
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

#ifndef RANGE_SPACE_STEP_STD_STEP_H
#define RANGE_SPACE_STEP_STD_STEP_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/AlgorithmStep.h"

namespace ReducedSpaceSQPPack {

///
/** Calculates the range space step by, solving for <tt>py = -inv(R)*c(con_decomp)</tt>,
 * then setting <tt>Ypy = Y * py</tt>.
 */
class RangeSpaceStepStd_Step
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

};	// end class RangeSpaceStepStd_Step

}	// end namespace ReducedSpaceSQPPack 

#endif	// RANGE_SPACE_STEP_STD_STEP_H
