// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateMultFree_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorClass.h"

namespace ReducedSpaceSQPPack {

MeritFunc_PenaltyParamUpdateMultFree_AddedStep::MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
		  const merit_func_ptr_t& merit_func, value_type small_mu
		, value_type mult_factor, value_type kkt_near_sol )
	: MeritFunc_PenaltyParamUpdateGuts_AddedStep(merit_func,small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateMultFree_AddedStep::min_mu(
	rSQPState& s, value_type* min_mu ) const
{
	using LinAlgPack::dot;
	using SparseLinAlgPack::dot;
	
	if ( s.Gf().updated_k(0) && s.nu().updated_k(0) && s.Ypy().updated_k(0) && s.c().updated_k(0) ) {
		// min_mu = abs((Gf_k+nu_k)'*Ypy_k) / norm(c_k,1)
		*min_mu = ::fabs( 	  dot( s.Gf().get_k(0)(), s.Ypy().get_k(0)() )
							+ dot( s.nu().get_k(0)(), s.Ypy().get_k(0)() )
						) / s.c().get_k(0).norm_1();
		return true;
	}
	return false;
}

void MeritFunc_PenaltyParamUpdateMultFree_AddedStep::print_min_mu_step(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if Gf_k, nu_k, Ypy_k and c_k are updated then\n"
		<< L << "    min_mu = abs((Gf_k+nu_k)'*Ypy_k) / norm(c_k,one)\n"
		<< L << "    update_mu = true\n"
		<< L << "else\n"
		<< L << "    update_mu = false\n"
		<< L << "endif\n"
		;
}

}	// end namespace ReducedSpaceSQPPack