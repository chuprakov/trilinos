// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_TangentialStepIP_Step.hpp
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

#ifndef NullSpaceStepIP_Step_H
#define NullSpaceStepIP_Step_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

///
/** Null Space Step for Interior Point algorithm
 *
 * This class calculates the step pz (and Zpz)
 *
 */

class TangentialStepIP_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
					 , poss_type assoc_step_poss);
		
		
		void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
						 , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
		//@}

	private:

	}; // end class TangentialStepIP_Step

} // end namespace MoochoPack

#endif // #if !defined NullSpaceStepIP_Step_H
