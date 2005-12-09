// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_SetDBoundsStd_AddedStep.hpp
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

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"

namespace MoochoPack {

///
/** Computes the bounds for the QP subproblem from the %NLP bounds.
 */
class SetDBoundsStd_AddedStep
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	///
	SetDBoundsStd_AddedStep();

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

private:
	IterationPack::CastIQMember<VectorMutable> dl_iq_;
	IterationPack::CastIQMember<VectorMutable> du_iq_;

}; // end class SetDBoundsStd_AddedStep

} // end namespace MoochoPack 

#endif // SET_D_BOUNDS_STD_ADDED_STEP_HH
