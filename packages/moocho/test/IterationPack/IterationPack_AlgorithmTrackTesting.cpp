// ///////////////////////////////////////////////////////////////
// AlgorithmTrackTesting.cpp

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <iomanip>

#include "../test/AlgorithmTrackTesting.h"
#include "../include/Algorithm.h"

namespace GeneralIterationPack {

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

}	// end namespace GeneralIterationPack 
