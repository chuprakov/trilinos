// ///////////////////////////////////////////////////////////////
// AlgorithmTrackTesting.h
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
