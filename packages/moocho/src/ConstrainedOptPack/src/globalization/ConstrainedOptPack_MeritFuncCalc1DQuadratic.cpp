// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalc1DQuadratic.cpp
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

#include "../include/MeritFuncCalc1DQuadratic.h"
#include "../include/MeritFuncCalc.h"
#include "LinAlgPack/include/VectorOp.h"

namespace ConstrainedOptimizationPack {

MeritFuncCalc1DQuadratic::MeritFuncCalc1DQuadratic( const MeritFuncCalc& phi, size_type p
		, const VectorSlice d[], VectorSlice* x )
	: phi_(phi), p_(p), x_(x)
{
	if( p < 1 || 2 < p )
		throw std::invalid_argument( "MeritFuncCalc1DQuadratic::MeritFuncCalc1DQuadratic(...) : "
			"Error, p must be in the range 1 <= p <= 2"								);
	for( size_type i = 0; i <= p; ++i )
		const_cast<VectorSlice&>(d_[i]).bind(const_cast<VectorSlice&>(d[i]));
}

value_type MeritFuncCalc1DQuadratic::operator()(value_type alpha) const
{
	using LinAlgPack::Vp_StV;
	*x_ = d_[0];
	value_type alpha_i = alpha;
	for( size_type i = 1; i <= p_; ++i, alpha_i *= alpha ) {
		Vp_StV( x_, alpha_i, d_[i] );
	}
	return phi_( *x_ );
}

value_type  MeritFuncCalc1DQuadratic::deriv() const
{
	return phi_.deriv();
}

void MeritFuncCalc1DQuadratic::print_merit_func( std::ostream& out
	, const std::string& L ) const
{
	out	<< L << "*** MeritFuncCalc1DQuadratic\n"
		<< L << "x = xo + alpha*d[1] + alpha^2*d[2]\n";
	phi_.print_merit_func( out, L );
}

}	// end namespace ConstrainedOptimizationPack
