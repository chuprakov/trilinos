// ///////////////////////////////////////////////////////////////
// AlgorithmTrack.h

#ifndef ALGORITHM_TRACK_H
#define ALGORITHM_TRACK_H

#include <iosfwd>

#include "GeneralIterationPackTypes.h"

namespace GeneralIterationPack {

///
/** Used to ouput iteration results and other information.
  *
  *
  *
  */
class AlgorithmTrack {
public:

	///
	AlgorithmTrack(std::ostream& journal_out)
		: journal_out_(&journal_out)
	{}

	///
	/** Output information about an iteration just completed.
	  *
	  * The default just does nothing.
	  */
	virtual void output_iteration(const Algorithm& algo) const {}

	///
	/** Output information about a just completed algorithm.
	  *
	  * The default just does nothing.
	  */
	virtual void output_final(const Algorithm& algo, EAlgoReturn algo_return) const {}

	///
	/** Return a reference to a std::ostream to be used to output debug information 
	  * and the like.
	  */
	virtual std::ostream& journal_out() const
	{	return *const_cast<AlgorithmTrack*>(this)->journal_out_; }

protected:
	std::ostream *journal_out_;

	// not defined and not to be called
	AlgorithmTrack();
	
};	// end class AlgorithmTrack

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_TRACK_H