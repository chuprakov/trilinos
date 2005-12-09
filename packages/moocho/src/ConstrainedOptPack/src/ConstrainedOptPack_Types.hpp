// ///////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_Types.hpp
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

#include "NLPInterfacePack_Types.hpp"
#include "NLPInterfacePack_NLP.hpp"

namespace ConstrainedOptPack {

#include "NLPInterfacePack_PublicTypes.ud"

/// Bounds type
enum EBounds { FREE, UPPER, LOWER, EQUALITY };

// concrete classes

class VariableBoundsTester;

// abstract classes

class MatrixSymAddDelUpdateableWithOpFactorized;
class MatrixIdentConcat;
class MeritFuncCalc1D;
class MeritFuncCalc;
class MeritFuncNLP;
class MeritFuncNLE;
class MeritFuncNLF;
class MeritFuncNLPDirecDeriv;
class MeritFuncPenaltyParam;
class MeritFuncPenaltyParams;
class DirectLineSearch_Strategy;

// concrete subclasses

class MeritFuncCalc1DQuadratic;
class MeritFuncCalcNLP;
class MeritFuncNLPL1;
class MeritFuncNLPModL1;
//class MeritFuncCalcNLE;
//class MeritFuncCalcNLF;
//class MatrixHessianSuperBasic;
//class MatrixHessianSuperBasicInitDiagonal;
//class MatrixSymPosDefInvCholFactor;
class MatrixSymPosDefLBFGS;
class MatrixSymAddDelBunchKaufman;
class MatrixSymHessianRelaxNonSing;
class MatrixIdentConcatStd;
class DirectLineSearchArmQuad_Strategy;
class DirectLineSearchArmQuad_StrategySetOptions;
class VarReductOrthogDenseStd_Strategy;

// decomposition classes

class DecompositionSystem;
class DecompositionSystemVarReduct;
class DecompositionSystemVarReductPerm;
class DecompositionSystemVarReductPermStd;
class DecompositionSystemVarReductImp;
class DecompositionSystemCoordinate;
class DecompositionSystemOrthogonal;
class DecompositionSystemTester;
class DecompositionSystemTesterSetOptions;

// Abstract QP solvers

class QPSolverRelaxed;
class QPSolverRelaxedTester;
class QPSolverRelaxedTesterSetOptions;

// Concrete QP solvers

//class QPSchur;
//class QPSolverRelaxedQPSchurRangeSpace;

}	// end namespace ConstrainedOptPack 

#endif // CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
