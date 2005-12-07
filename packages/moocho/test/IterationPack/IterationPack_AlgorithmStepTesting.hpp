// ///////////////////////////////////////////////////////////////
// IterationPack_AlgorithmStepTesting.hpp
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

#ifndef ALGORITHM_STEP_TESTING_H
#define ALGORITHM_STEP_TESTING_H

#include "IterationPack_AlgorithmStep.hpp"

namespace IterationPack {

///
/** Testing class
  */
class AlgorithmStepTesting : public AlgorithmStep {
public:

	///
	bool do_step(
		Algorithm& algo, poss_type step_poss, EDoStepType type
		,poss_type assoc_step_poss
		);
	///
	void initialize_step(
		Algorithm& algo, poss_type step_poss, EDoStepType type
		,poss_type assoc_step_poss
		);
	///
	void inform_updated(
		Algorithm& algo, poss_type step_poss, EDoStepType type
		,poss_type assoc_step_poss
		);
	///
	void finalize_step(
		Algorithm& algo, poss_type step_poss, EDoStepType type
		,poss_type assoc_step_poss
		);
	///
	void print_step(
		const Algorithm& algo, poss_type step_poss, EDoStepType type
		,poss_type assoc_step_poss ,std::ostream& out, const std::string& leading_str
		) const;

private:
	
	///
	void print_step_poss(
		const Algorithm& algo, poss_type step_poss, EDoStepType type
		,poss_type assoc_step_poss
		) const;

};	// end class AlgorithmStepTesting

}	// end namespace IterationPack 

#endif // ALGORITHM_STEP_TESTING_H
