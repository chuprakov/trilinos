// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_MeritFunc_DummyUpdate_Step.hpp
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

#ifndef MERIT_FUNC_DUMMY_UPDATE_STEP_H
#define MERIT_FUNC_DUMMY_UPDATE_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

///
/** Simply updates <tt>merit_func_nlp_k = merit_func_nlp_km1</tt>
  */
class MeritFunc_DummyUpdate_Step
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

};	// end class MeritFunc_DummyUpdate_Step

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_DUMMY_UPDATE_STEP_H
