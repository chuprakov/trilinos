// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLPL1.cpp
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

#include "../include/MeritFuncNLPL1.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"

namespace ConstrainedOptimizationPack {

MeritFuncNLPL1::MeritFuncNLPL1()
	: deriv_(0.0), mu_(0.0)
{}

// Overridden from MeritFuncNLP

value_type MeritFuncNLPL1::value(value_type f, const VectorSlice& c) const
{
	using LinAlgPack::norm_1;
	return f + mu_ * norm_1(c);
}

value_type MeritFuncNLPL1::deriv() const
{
	return deriv_;
}

void MeritFuncNLPL1::print_merit_func(std::ostream& out
	, const std::string& L ) const
{
	out
		<< L << "*** Define L1 merit funciton (assumes Gc_k'*d_k + c_k = 0):\n"
		<< L << "phi(f,c) = f + mu_k * norm(c,1)\n"
		<< L << "Dphi(x_k,d_k) = Gf_k' * d_k - mu * norm(c_k,1)\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPL1::calc_deriv( const VectorSlice& Gf_k, const VectorSlice& c_k
	, const VectorSlice& d_k )
{
	using LinAlgPack::dot; using LinAlgPack::norm_1;
	return deriv_ = dot( Gf_k, d_k ) - mu_ * norm_1( c_k );
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPL1::mu(value_type mu)
{
	mu_ = mu;
}

value_type MeritFuncNLPL1::mu() const
{
	return mu_;
}

}	// end namespace ConstrainedOptimizationPack 
