// /////////////////////////////////////////////////////////////////////////////
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

#include "ConstrainedOptimizationPack/include/MeritFuncNLPL1.h"
#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "ThrowException.h"

namespace ConstrainedOptimizationPack {

MeritFuncNLPL1::MeritFuncNLPL1()
	: deriv_(0.0), mu_(0.0)
{}

// Overridden from MeritFuncNLP

value_type MeritFuncNLPL1::value(
	value_type             f
	,const VectorWithOp    *c
	,const VectorWithOp    *h
	,const VectorWithOp    *hl
	,const VectorWithOp    *hu
	) const
{
	THROW_EXCEPTION(
		h || hl || hu, std::logic_error
		,"MeritFuncNLPL1::value(...) : Error! general inequalities are not supported yet" );
	return f + ( c ? mu_ * c->norm_1() : 0.0);
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

value_type MeritFuncNLPL1::calc_deriv(
	const VectorWithOp    &Gf_k
	,const VectorWithOp   *c_k
	,const VectorWithOp   *h_k
	,const VectorWithOp   *hl
	,const VectorWithOp   *hu
	,const VectorWithOp   &d_k
	)
{
	using AbstractLinAlgPack::dot;
	THROW_EXCEPTION(
		h_k || hl || hu, std::logic_error
		,"MeritFuncNLPL1::value(...) : Error! general inequalities are not supported yet" );
	return deriv_ = dot( Gf_k, d_k ) - ( c_k ? mu_ * c_k->norm_1() : 0.0 );
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
