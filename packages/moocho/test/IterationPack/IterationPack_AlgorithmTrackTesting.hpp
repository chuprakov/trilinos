// ///////////////////////////////////////////////////////////////
// AlgorithmTrackTesting.h

#ifndef ALGORITHM_TRACK_TESTING_H
#define ALGORITHM_TRACK_TESTING_H

#include "../include/AlgorithmTrack.h"

namespace GeneralIterationPack {

///
/** Testing class.
  */
class AlgorithmTrackTesting : public AlgorithmTrack {
public:

	AlgorithmTrackTesting(std::ostream& out) : AlgorithmTrack(out)
	{}

	// Overriden
	
	///
	void output_iteration(const Algorithm& algo) const;

	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

};	// end class AlgorithmTrackTesting

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_TRACK_TESTING_H