// ////////////////////////////////////////////////////////////////////////
// MeritFuncNLPModL1.cpp
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

#include "ConstrainedOptimizationPack/src/MeritFuncNLPModL1.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorStdOps.hpp"
#include "ThrowException.hpp"

namespace ConstrainedOptimizationPack {

MeritFuncNLPModL1::MeritFuncNLPModL1()
	: deriv_(0.0)
{}

// Overridden from MeritFuncNLP

value_type MeritFuncNLPModL1::value(
	value_type             f
	,const VectorWithOp    *c
	,const VectorWithOp    *h
	,const VectorWithOp    *hl
	,const VectorWithOp    *hu
	) const
{
	THROW_EXCEPTION(
		h || hl || hu, std::logic_error
		,"MeritFuncNLPModL1::value(...) : Error! general inequalities are not supported!" );
/*
	using LinAlgPack::norm_1;
	return f + local_constr_term( mu_, c, "calc_deriv" );
*/
	assert(0); // ToDo: Write a reduction operator for the above operation
	return 0.0;
}

value_type MeritFuncNLPModL1::deriv() const
{
	return deriv_;
}

void MeritFuncNLPModL1::print_merit_func(
	std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Define a modified L1 merit funciton that uses different\n"
		<< L << "*** penalty parameters for each constriant.\n"
		<< L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
		<< L << "phi(f,c) = f + sum( mu(j) * abs(c(j)), j = 1,...,m )\n"
		<< L << "Dphi(x_k,d_k) = Gf_k' * d_k - sum( mu(j) * abs(c(j)), j = 1,...,m )\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPModL1::calc_deriv(
	const VectorWithOp    &Gf_k
	,const VectorWithOp   *c_k
	,const VectorWithOp   *h_k
	,const VectorWithOp   *hl
	,const VectorWithOp   *hu
	,const VectorWithOp   &d_k
	)
{
	THROW_EXCEPTION(
		h_k || hl || hu, std::logic_error
		,"MeritFuncNLPModL1::value(...) : Error! general inequalities are not supported!" );
/*
	using LinAlgPack::dot; using LinAlgPack::norm_1;
	return deriv_ = dot( Gf_k, d_k ) - local_constr_term( mu_, c_k, "calc_deriv" );
*/
	assert(0); // ToDo: Write a reduction operator for the above operation
	return 0.0;
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPModL1::set_space_c( const VectorSpace::space_ptr_t& space_c )
{
	mu_  = space_c->create_member();
	*mu_ = 0.0;
}

VectorWithOpMutable& MeritFuncNLPModL1::set_mu()
{
	return *mu_;
}

const VectorWithOp& MeritFuncNLPModL1::get_mu() const
{
	return *mu_;
}

}	// end namespace ConstrainedOptimizationPack

/* ToDo: Write a reduction operator for the following!

namespace {

value_type local_constr_term( const Vector& mu, const VectorSlice& c
	, const char func_name[] )
{
	if( mu.size() != c.size() ) {
		std::ostringstream omsg;
		omsg
			<< "MeritFuncNLPModL1::" << func_name << "(...) : "
			<< "Error, the sizes mu.size() == " << mu.size()
			<< " != c.size() == " << c.size();
		throw ConstrainedOptimizationPack::MeritFuncNLP::InvalidInitialization(omsg.str());
	}
	value_type r = 0.0;
	Vector::const_iterator
		mu_itr = mu.begin();
	VectorSlice::const_iterator
		c_itr = c.begin();
	while( mu_itr != mu.end() )
		r += *mu_itr++ * ::fabs( *c_itr++ );
	return r;
}

}	// end namespace

*/
