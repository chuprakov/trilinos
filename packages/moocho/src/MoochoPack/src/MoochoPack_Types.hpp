// ///////////////////////////////////////////////////////////////////////
// ReducedSpaceSQPPackTypes.h
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

#ifndef REDUCED_SPACE_SQP_PACK_TYPES_H
#define REDUCED_SPACE_SQP_PACK_TYPES_H

#include "ConstrainedOptimizationPack/include/ConstrainedOptimizationPackTypes.h"
#include "GeneralIterationPack/include/GeneralIterationPackTypes.h"

namespace ReducedSpaceSQPPack {

// using types from ConstrainedOptimizationPack
#include "ConstrainedOptimizationPack/include/ConstrainedOptimizationPackPublicTypes.ud"

// using types from GeneralIterationPack
#include "GeneralIterationPack/include/GeneralIterationPackPublicTypes.ud"

// enum for rSQP output.
enum EJournalOutputLevel {
	PRINT_NOTHING = 0,
	PRINT_BASIC_ALGORITHM_INFO = 1,
	PRINT_ALGORITHM_STEPS = 2,
	PRINT_ACTIVE_SET = 3,
	PRINT_VECTORS = 4,
	PRINT_ITERATION_QUANTITIES = 5
};

// public interface classes

class rSQPState;
class rSQPTrack;
class rSQPSolverClientInterface;
class rSQPAlgoClientInterface;
class rSQPAlgo_Config;
class IterQuantMatrixWithOpCreator;

//

class rSQPAlgo;

}	// end namespace ReducedSpaceSQPPack 

#endif // REDUCED_SPACE_SQP_PACK_TYPES_H
