// ///////////////////////////////////////////////////////////////
// AlgorithmTrackTesting.cpp
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

#include <iomanip>

#include "../test/AlgorithmTrackTesting.hpp"
#include "IterationPack/src/Algorithm.hpp"

namespace IterationPack {

void AlgorithmTrackTesting::output_iteration(const Algorithm& algo) const {
	journal_out()
			<< "\ntrack.output_iteration(algo) called for the iteration k = "
			<< algo.state().k() << std::endl;
}

void AlgorithmTrackTesting::output_final(const Algorithm& algo, EAlgoReturn algo_return) const {
	char algo_return_name[3][20] = { "TERMINATE_TRUE" ,"TERMINATE_FALSE" ,"MAX_ITER_EXCEEDED" };
	journal_out()
			<< "\ntrack.output_final(algo,algo_return) called for the iteration k = "
			<< algo.state().k() << " and algo_return = " << algo_return_name[algo_return]
			<< std::endl;
}

}	// end namespace IterationPack 
