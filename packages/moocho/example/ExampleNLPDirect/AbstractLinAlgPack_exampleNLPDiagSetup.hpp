// ////////////////////////////////////////////////////////////////////
// exampleNLPDiagSetup.hpp

#ifndef ALAP_EXPL_NLP_DIAG_SETUP_HPP
#define ALAP_EXPL_NLP_DIAG_SETUP_HPP

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ExampleVectorLib/src/MPIDenseVector.hpp"
#include "AbstractLinAlgPack/src/serial/interfaces/VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack/src/abstract/tsfl/VectorSpaceTSFL.hpp"
#include "TSFCore/src/Core/VectorSpaceSerialDecl.hpp"
#include "CommandLineProcessor.hpp"

namespace AbstractLinAlgPack {

///
/** Create a vector space given the input arguments argc, argv[] and an MPI communicator
 *
 */
int exampleNLPDiagSetup(
	int argc, char* argv[], MPI_Comm comm
	,MemMngPack::ref_count_ptr<const VectorSpace> *vec_space
	,size_type *n, value_type *xo, bool *has_bounds, bool *dep_bounded
	);

} // namespace AbstractLinAlgPack

#endif // ALAP_EXPL_NLP_DIAG_SETUP_HPP
