// ///////////////////////////////////////////////////////////////////////
// ConstrainedOptimizationPackTypes.h
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

#ifndef CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
#define CONSTRAINED_OPTIMIZATION_PACK_TYPES_H

#include "ConstrainedOptimizationPackDebugAcronyms.h"
#include "SparseLinAlgPack/include/SparseLinAlgPackTypes.h"
//#include "SparseSolverPack/include/SparseSolverPackTypes.h"
#include "NLPInterfacePack/include/NLP.h"

namespace ConstrainedOptimizationPack {

#include "SparseLinAlgPack/include/SparseLinAlgPackPublicTypes.ud"
//#include "SparseSolverPack/include/SparseSolverPackPublicTypes.ud"
#include "NLPInterfacePack/include/NLPInterfacePackPublicTypes.ud"

/// Bounds type
enum EBounds { FREE, UPPER, LOWER, EQUALITY };

// concrete classes

/*
class VectorWithNorms;
class DenseIdentVertConcatMatrix;
class IdentZeroVertConcatMatrix;
*/

// abstract classes

/*
class MatrixSymSecantUpdateable;
class MatrixSymAddDelUpdateable;
class MatrixSymAddDelUpdateableWithOpFactorized;
class MeritFuncCalc1D;
class MeritFuncCalc;
class MeritFuncNLP;
class MeritFuncNLE;
class MeritFuncNLF;
class MeritFuncNLPDirecDeriv;
class MeritFuncPenaltyParam;
class MeritFuncPenaltyParams;
class DirectLineSearch_Strategy;
class ZVarReductMatrix;
*/

// concrete subclasses

/*
class MeritFuncCalc1DQuadratic;
class MeritFuncCalcNLP;
class MeritFuncCalcNLE;
class MeritFuncCalcNLF;
class MatrixHessianSuperBasic;
class MatrixHessianSuperBasicInitDiagonal;
class MatrixSymPosDefInvCholFactor;
class MatrixSymPosDefCholFactor;
class MatrixSymPosDefLBFGS;
class MatrixSymAddDelBunchKaufman;
class DenseIdentVertConcatMatrixSubclass;
class ZAdjointFactMatrixSubclass;
class IdentZeroVertConcatMatrixSubclass;
*/

// decomposition classes

/*
class DecompositionSystem;
class DecompositionSystemVarReduct;
class DecompositionSystemVarReductImpNode;
class DecompositionSystemCoordinate;
class DecompositionSystemCoordinateDirect;
class DecompositionSystemCoordinateAdjoint;
*/

// Abstract QP solvers

//class QPSolverWithBounds;
//class QPSolverRelaxed;

// Concrete QP solvers

//class QPSCPD;
//class QPSchur;
//class QPSolverRelaxedQPSchurRangeSpace;

}	// end namespace ConstrainedOptimizationPack 

#endif // CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
