// ////////////////////////////////////////////////////////////////////////////
// SetDBoundsStd_AddedStep.h
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

#ifndef SET_D_BOUNDS_STD_ADDED_STEP_HH
#define SET_D_BOUNDS_STD_ADDED_STEP_HH

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/src/AlgorithmStep.h"
#include "GeneralIterationPack/src/CastIQMember.h"
#include "ReducedSpaceSQPPack/src/std/d_bounds_iter_quant.h"

namespace ReducedSpaceSQPPack {

///
/** Computes the bounds for the QP subproblem from the %NLP bounds.
 */
class SetDBoundsStd_AddedStep
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
{
public:

	///
	SetDBoundsStd_AddedStep();

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

private:
	GeneralIterationPack::CastIQMember<VectorWithOpMutable> dl_iq_;
	GeneralIterationPack::CastIQMember<VectorWithOpMutable> du_iq_;

}; // end class SetDBoundsStd_AddedStep

} // end namespace ReducedSpaceSQPPack 

#endif // SET_D_BOUNDS_STD_ADDED_STEP_HH
