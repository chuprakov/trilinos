// ///////////////////////////////////////////////////////////////////////
// ReducedSpaceSQPPackTypes.hpp
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

#include "ConstrainedOptimizationPack/src/ConstrainedOptimizationPackTypes.hpp"
#include "IterationPack/src/IterationPackTypes.hpp"

namespace ReducedSpaceSQPPack {

// using types from ConstrainedOptimizationPack
#include "ConstrainedOptimizationPack/src/ConstrainedOptimizationPackPublicTypes.ud"

// using types from IterationPack
#include "IterationPack/src/IterationPackPublicTypes.ud"

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
class rSQPSolverClientInterface;
class rSQPAlgoClientInterface;
class rSQPAlgo_Config;

//

class rSQPAlgo;
typedef IterationPack::AlgorithmStep             rSQPStep;
typedef IterationPack::AlgorithmTracker            rSQPTrack;
typedef IterationPack::AlgorithmTrackerComposite   rSQPTrackComposite;

}	// end namespace ReducedSpaceSQPPack 

#endif // REDUCED_SPACE_SQP_PACK_TYPES_H
