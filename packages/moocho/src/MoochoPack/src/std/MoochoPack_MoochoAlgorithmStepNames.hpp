// //////////////////////////////////////////////////////////////
// MoochoPack_MoochoAlgorithmStepNames.hpp
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

#ifndef RSQP_ALGORITHM_STEP_NAMES_H
#define RSQP_ALGORITHM_STEP_NAMES_H

#include <string>

#include "../MoochoPack_Types.hpp"

namespace MoochoPack {

/** @name Names for MOOCHOsteps */
//@{

const std::string EvalNewPoint_name                 = "EvalNewPoint";
const std::string ReducedGradient_name              = "ReducedGradient";
const std::string ReducedHessian_name               = "ReducedHessian";
const std::string QuasiNormalStep_name              = "QuasiNormalStep";
const std::string TangentialStep_name               = "TangentialStep";
const std::string SearchDirec_name                  = "SearchDirec";
const std::string LineSearch_name                   = "LineSearch";
const std::string CheckConvergence_name             = "CheckConvergence";

const std::string CalcLambdaIndep_name              = "CalcLambdaIndep";
const std::string CalcReducedGradLagrangian_name    = "CalcReducedGradLagrangian";
const std::string CheckSkipBFGSUpdate_name          = "CheckSkipBFGSUpdate";
const std::string CalcDFromYPYZPZ_name				= "CalcDFromYPYZPZ";

//@}
}	// end namespace MoochoPack

#endif // RSQP_ALGORITHM_STEP_NAMES_H
