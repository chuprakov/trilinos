// ///////////////////////////////////////////////////////////////
// AlgorithmStep.h

#ifndef ALGORITHM_STEP_H
#define ALGORITHM_STEP_H

#include <iosfwd>

#include "GeneralIterationPackTypes.h"

namespace GeneralIterationPack {

///
/** Base type for all objects that perform steps in an #Algorithm#.
  *
  *
  *
  */
class AlgorithmStep {
public:

	///
	typedef size_t poss_type;

	///
	/** Called by #Algorithm# to perform a main, pre or post step at step_poss and assoc_step_poss.
	  *
	  */
	virtual bool do_step( Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss ) = 0;

	///
	/** Called by #Algorithm# to inform when a runtime configuration change
	  * is finihed.
	  *
	  * The default does nothing.
	  */
	virtual void inform_updated(Algorithm& algo)
	{}

	///
	/** Called by #Algorithm# to print out what this step does in Matlab like format.
	  *
	  * The default does nothing.
	  */
	virtual void print_step( const Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const
	{}

};	// end class AlgorithmStep

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_STEP_H