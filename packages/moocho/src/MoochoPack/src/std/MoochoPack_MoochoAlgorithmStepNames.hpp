// //////////////////////////////////////////////////////////////
// rSQPAlgorithmStepNames.h

#ifndef RSQP_ALGORITHM_STEP_NAMES_H
#define RSQP_ALGORITHM_STEP_NAMES_H

#include <string>

#include "../ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

/** @name Names for rSQP++ steps */
//@{

const std::string EvalNewPoint_name					= "EvalNewPoint";
const std::string ReducedGradient_name				= "ReducedGradient";
const std::string ReducedHessian_name				= "ReducedHessian";
const std::string DepDirec_name						= "DepDirec";
const std::string IndepDirec_name					= "IndepDirec";
const std::string SearchDirec_name					= "SearchDirec";
const std::string LineSearch_name					= "LineSearch";
const std::string CheckConvergence_name				= "CheckConvergence";

const std::string CalcLambdaIndep_name				= "CalcLambdaIndep";
const std::string CalcReducedGradLagrangian_name	= "CalcReducedGradLagrangian";
const std::string CheckSkipBFGSUpdate_name			= "CheckSkipBFGSUpdate";

//@}
}	// end namespace ReducedSpaceSQPPack

#endif // RSQP_ALGORITHM_STEP_NAMES_H
