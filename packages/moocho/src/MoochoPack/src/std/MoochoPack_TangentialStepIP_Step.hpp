// ////////////////////////////////////////////////////////////////////////////
// NullSpaceStepIP_Step.hpp
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

#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackTypes.hpp"
#include "GeneralIterationPack/src/AlgorithmStep.hpp"

namespace ReducedSpaceSQPPack {

///
/** Null Space Step for Interior Point algorithm
 *
 * This class calculates the step pz (and Zpz)
 *
 */

class NullSpaceStepIP_Step
	: public GeneralIterationPack::AlgorithmStep // doxygen needs full path
	{
	public:

		/** @name Overridden from AlgorithmStep */
		//@{
		///
		bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
					 , poss_type assoc_step_poss);
		
		
		void print_step( const GeneralIterationPack::Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
						 , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
		//@}

	private:

	}; // end class NullSpaceStepIP_Step

} // end namespace ReducedSpaceSQPPack

#endif // #if !defined NullSpaceStepIP_Step_H
