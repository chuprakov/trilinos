// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_CalcD_vStep_Step.hpp
//
// Copyright (C) 2001
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

#ifndef CALC_D_v_STEP_STEP_H
#define CALC_D_v_STEP_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

///
/** Calculates <tt>dvl_k = mu*invXl_k*e - vl_k - invXl_k*Vl_k*d_k</tt>
 *   and       <tt>dvu_k = mu*invXu_k*e - vu_k + invXu_k*Vu_k*d_k</tt>
 */
class CalcD_vStep_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

};	// end class CalcD_vStep_Step

}	// end namespace MoochoPack 

#endif	// CALC_D_v_STEP_STEP_H
