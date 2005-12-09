// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalcNLE.cpp
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

#include "ConstrainedOptPack_MeritFuncCalcNLE.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

namespace ConstrainedOptPack {

MeritFuncCalcNLE::MeritFuncCalcNLE( const MeritFuncNLE* phi, const NLP* nlp )
	: phi_(phi), nlp_(nlp)
{}

value_type MeritFuncCalcNLE::operator()(const Vector& x) const {
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

}	// end namespace ConstrainedOptPack
