// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalcNLE.cpp

#include "../include/MeritFuncCalcNLE.h"
#include "LinAlgPack/include/VectorClass.h"

namespace ConstrainedOptimizationPack {

MeritFuncCalcNLE::MeritFuncCalcNLE( const MeritFuncNLE* phi, const NLP* nlp )
	: phi_(phi), nlp_(nlp)
{}

value_type MeritFuncCalcNLE::operator()(const VectorSlice& x) const {
	nlp().calc_c(x);
	return phi().value( nlp().c() );
}

value_type MeritFuncCalcNLE::deriv() const {
	return phi().deriv();
}

void MeritFuncCalcNLE::print_merit_func(std::ostream& out
	, const std::string& L) const
{
	out	<< L << "*** MeritFuncCalcNLE\n"
		<< L << "c = c(x)\n";
	phi().print_merit_func(out,L);
}

}	// end namespace ConstrainedOptimizationPack
