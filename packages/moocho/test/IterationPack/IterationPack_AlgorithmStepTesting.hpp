// ///////////////////////////////////////////////////////////////
// AlgorithmStepTesting.h

#ifndef ALGORITHM_STEP_TESTING_H
#define ALGORITHM_STEP_TESTING_H

#include "../include/AlgorithmStep.h"

namespace GeneralIterationPack {

///
/** Testing class
  */
class AlgorithmStepTesting : public AlgorithmStep {
public:

	///
	bool do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss);

	///
	void inform_updated(Algorithm& algo);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss ,std::ostream& out, const std::string& leading_str ) const;

};	// end class AlgorithmStepTesting

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_STEP_TESTING_H