// ///////////////////////////////////////////////////////////////////////
// ConstrainedOptimizationPackTypes.h

#ifndef CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
#define CONSTRAINED_OPTIMIZATION_PACK_TYPES_H

#include "ConstrainedOptimizationPackDebugAcronyms.h"
#include "SparseSolverPack/include/SparseSolverPackTypes.h"
#include "NLPInterfacePack/include/NLP.h"

namespace ConstrainedOptimizationPack {

// types from SparseSolverPack
#include "SparseSolverPack/include/SparseSolverPackPublicTypes.ud"

// types from NLPInterfacePack
#include "NLPInterfacePack/include/NLPInterfacePackPublicTypes.ud"

/// Bounds type
enum EBounds { FREE, UPPER, LOWER, EQUALITY };

// concrete classes

class VectorWithNorms;
class DenseIdentVertConcatMatrix;
class IdentZeroVertConcatMatrix;

// abstract classes

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

// concrete subclasses

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

// decomposition classes

class DecompositionSystem;
class DecompositionSystemVarReduct;
class DecompositionSystemVarReductImpNode;
class DecompositionSystemCoordinate;
class DecompositionSystemCoordinateDirect;
class DecompositionSystemCoordinateAdjoint;

// Abstract QP solvers
class QPSolverWithBounds;
class QPSolverRelaxed;

// Concrete QP solvers
class QPSCPD;
class QPSchur;
class QPSolverRelaxedQPSchurRangeSpace;

}	// end namespace ConstrainedOptimizationPack 

#endif // CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
