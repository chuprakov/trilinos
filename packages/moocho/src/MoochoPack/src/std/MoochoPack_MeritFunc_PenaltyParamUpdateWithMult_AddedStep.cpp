// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateWithMult_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateWithMult_AddedStep.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"

namespace ReducedSpaceSQPPack {

MeritFunc_PenaltyParamUpdateWithMult_AddedStep::MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
		  const merit_func_ptr_t& merit_func, value_type small_mu
		, value_type mult_factor, value_type kkt_near_sol )
	: MeritFunc_PenaltyParamUpdateGuts_AddedStep(merit_func,small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateWithMult_AddedStep::min_mu(
	rSQPState& s, value_type* min_mu ) const
{
	if ( s.lambda().updated_k(0) ) {
		*min_mu = s.lambda().get_k(0).norm_inf();
		return true;
	}
	return false;
}

void MeritFunc_PenaltyParamUpdateWithMult_AddedStep::print_min_mu_step(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if lambda_k is updated then\n"
		<< L << "    min_mu = norm( lambda_k, inf )\n"
		<< L << "    update_mu = true\n"
		<< L << "else\n"
		<< L << "    update_mu = false\n"
		<< L << "endif\n"
		;
}

}	// end namespace ReducedSpaceSQPPack
