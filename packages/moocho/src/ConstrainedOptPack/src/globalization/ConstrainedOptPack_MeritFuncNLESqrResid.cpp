// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLESqrResid.cpp

#include "../include/MeritFuncNLESqrResid.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"

namespace ConstrainedOptimizationPack {

MeritFuncNLESqrResid::MeritFuncNLESqrResid()
	: deriv_(0.0)
{}

value_type MeritFuncNLESqrResid::calc_deriv( const VectorSlice& c_k )
{
	using LinAlgPack::dot;
	return deriv_ = - dot(c_k,c_k);
}

// Overridden from MeritFuncNLP

value_type MeritFuncNLESqrResid::value(const VectorSlice& c) const
{
	using LinAlgPack::dot;
	return 0.5 * dot(c,c);
}

value_type MeritFuncNLESqrResid::deriv() const
{
	return deriv_;
}

void MeritFuncNLESqrResid::print_merit_func(std::ostream& out
	, const std::string& L ) const
{
	out
		<< L << "*** Define a square of constraint residuals merit funciton\n"
		<< L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
		<< L << "phi(c) = 1/2 * dot(c,c)\n"
		<< L << "Dphi(x_k,d_k) = - dot(c_k,c_k)\n";
}

}	// end namespace ConstrainedOptimizationPack 
