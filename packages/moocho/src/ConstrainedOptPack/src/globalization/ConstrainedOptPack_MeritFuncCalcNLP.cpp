// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalcNLP.cpp
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

#include "MeritFuncCalcNLP.hpp"
#include "LinAlgPack/src/VectorClass.hpp"

namespace ConstrainedOptimizationPack {

MeritFuncCalcNLP::MeritFuncCalcNLP( const MeritFuncNLP* phi, const NLP* nlp )
	: phi_(phi), nlp_(nlp)
{}

value_type MeritFuncCalcNLP::operator()(const VectorWithOp& x) const
{
	const size_type
		m  = nlp().m(),
		mI = nlp().mI();
	nlp().calc_f(x);
	if(m)  nlp().calc_c(x,false);
	if(mI) nlp().calc_h(x,false);
	return phi().value(
		nlp().f()
		,m  ? &nlp().c()  : NULL
		,mI ? &nlp().h()  : NULL
		,mI ? &nlp().hl() : NULL
		,mI ? &nlp().hu() : NULL
		);
}

value_type MeritFuncCalcNLP::deriv() const {
	return phi().deriv();
}

void MeritFuncCalcNLP::print_merit_func(
	std::ostream& out, const std::string& L
	) const
{
	out	<< L << "*** MeritFuncCalcNLP\n"
		<< L << "f = f(x), c = c(x)\n";
	phi().print_merit_func(out,L);
}

}	// end namespace ConstrainedOptimizationPack
