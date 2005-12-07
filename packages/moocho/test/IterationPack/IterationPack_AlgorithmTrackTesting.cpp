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

#include "IterationPack_AlgorithmTrackTesting.hpp"
#include "IterationPack_Algorithm.hpp"

namespace IterationPack {

void AlgorithmTrackTesting::output_iteration(const Algorithm& algo) const {
	std::ostream &o = journal_out();
	o
		<< "\ntrack.output_iteration(algo) called for the iteration k = "
		<< algo.state().k() << std::endl;
	o
		<< "\nStep times for the current iteration (in seconds):\n";
	const int n = algo.num_steps();
	std::vector<double> step_times(n+1);
	algo.get_step_times_k(0,&step_times[0]);
	o << "  step_id:time = ";
	for( int l = 0; l < n+1; ++l )
		o << " " << (l+1) << ":" << step_times[l];
	o << "\n  total time = " << step_times[n] << std::endl;
	
}

void AlgorithmTrackTesting::output_final(const Algorithm& algo, EAlgoReturn algo_return) const {
	char algo_return_name[6][50] =
		{
			"TERMINATE_TRUE"
			,"TERMINATE_FALSE"
			,"MAX_ITER_EXCEEDED"
			,"MAX_RUN_TIME_EXCEEDED"
			,"INTERRUPTED_TERMINATE_TRUE"
			,"INTERRUPTED_TERMINATE_FALSE"
		};
	std::ostream &o = journal_out();
	o << "\ntrack.output_final(algo,algo_return) called for the iteration k = "
	  << algo.state().k() << " and algo_return = " << algo_return_name[algo_return]
	  << std::endl;
	o << "Timing (in seconds) statistics for step 0 : ";
	double total, average, min, max, percent;
	algo.get_final_step_stats(0,&total,&average,&min,&max,&percent);
	o << "total = " << total << ", average = " << average << ", min = " << min
	  << ", max = " << max << ", percent = " << percent << std::endl;
}

}	// end namespace IterationPack 
