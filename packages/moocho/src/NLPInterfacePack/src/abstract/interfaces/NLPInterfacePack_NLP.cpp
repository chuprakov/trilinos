// ///////////////////////////////////////////////////////////////////////
// NLP.cpp
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

#include <limits>

#include "NLPInterfacePack/include/NLP.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "ThrowException.h"

namespace {
const char name_f[] = "f";
const char name_c[] = "c";
const char name_h[] = "h";
NLPInterfacePack::NLP::options_ptr_t  null_options = MemMngPack::null;
} // end namespace

// static

NLPInterfacePack::value_type NLPInterfacePack::NLP::infinite_bound()
{
	return std::numeric_limits<value_type>::max();
	//	return 1e+50;
}

// constructors

NLPInterfacePack::NLP::NLP()
{}

// destructor

NLPInterfacePack::NLP::~NLP()
{}

void NLPInterfacePack::NLP::set_options( const options_ptr_t& options )
{}

const NLPInterfacePack::NLP::options_ptr_t&
NLPInterfacePack::NLP::get_options() const
{
	return null_options;
}

void NLPInterfacePack::NLP::initialize(bool test_setup)
{
	num_f_evals_ = num_c_evals_ = num_h_evals_ = 0;
}

// dimensionality

NLPInterfacePack::size_type
NLPInterfacePack::NLP::n() const
{
	return this->space_x()->dim();
}

NLPInterfacePack::size_type 
NLPInterfacePack::NLP::m() const
{
	VectorSpace::space_ptr_t spc = this->space_c();
	return spc.get() ? spc->dim() : 0;
}


NLPInterfacePack::size_type
NLPInterfacePack::NLP::mI() const
{
	VectorSpace::space_ptr_t spc = this->space_h();
	return spc.get() ? spc->dim() : 0;
}

// initial guess

void NLPInterfacePack::NLP::get_init_lagrange_mult(
	VectorWithOpMutable*   lambda
	,VectorWithOpMutable*  lambdaI
	,VectorWithOpMutable*  nu
	) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( lambda  && this->m()  == 0,            std::logic_error, "" );
	THROW_EXCEPTION( lambdaI && this->mI() == 0,            std::logic_error, "" );
	THROW_EXCEPTION( nu      && this->num_bounded_x() == 0, std::logic_error, "" );
#endif
	if(lambda) {
#ifdef _DEBUG
		THROW_EXCEPTION( !this->space_c()->is_compatible(lambda->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
		*lambda = 0.0;
	}
	if(lambdaI) {
#ifdef _DEBUG
		THROW_EXCEPTION( !this->space_h()->is_compatible(lambdaI->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
		*lambdaI = 0.0;
	}
	if(nu) {
#ifdef _DEBUG
		THROW_EXCEPTION( !this->space_x()->is_compatible(nu->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
		*nu = 0.0;
	}
}

// <<std comp>> members for f

void NLPInterfacePack::NLP::set_f(value_type* f)
{
	first_order_info_.f = f;
}

NLPInterfacePack::value_type* NLPInterfacePack::NLP::get_f()
{
	return StandardCompositionRelationshipsPack::get_role_name(first_order_info_.f, false, name_f);
}

NLPInterfacePack::value_type& NLPInterfacePack::NLP::f()
{
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.f, false, name_f);
}

const NLPInterfacePack::value_type& NLPInterfacePack::NLP::f() const
{
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.f, false, name_f);
}

// <<std comp>> members for c

void NLPInterfacePack::NLP::set_c(VectorWithOpMutable* c)
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
	THROW_EXCEPTION( c && !this->space_c()->is_compatible(c->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
	first_order_info_.c = c;
}

NLPInterfacePack::VectorWithOpMutable* NLPInterfacePack::NLP::get_c()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(first_order_info_.c, false, name_c);
}

NLPInterfacePack::VectorWithOpMutable& NLPInterfacePack::NLP::c()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.c, false, name_c);
}

const NLPInterfacePack::VectorWithOp& NLPInterfacePack::NLP::c() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.c, false, name_c);
}

// <<std comp>> members for h

void NLPInterfacePack::NLP::set_h(VectorWithOpMutable* h)
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
	THROW_EXCEPTION( h && !this->space_h()->is_compatible(h->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
	first_order_info_.h = h;
}

NLPInterfacePack::VectorWithOpMutable* NLPInterfacePack::NLP::get_h()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(first_order_info_.h, false, name_h);
}

NLPInterfacePack::VectorWithOpMutable& NLPInterfacePack::NLP::h()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.h, false, name_h);
}

const NLPInterfacePack::VectorWithOp& NLPInterfacePack::NLP::h() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.h, false, name_h);
}

// calculations

void NLPInterfacePack::NLP::set_multi_calc(bool multi_calc) const
{}

bool NLPInterfacePack::NLP::multi_calc() const
{
	return false;
}

void NLPInterfacePack::NLP::calc_f(const VectorWithOp& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_.f, "NLP::calc_f()", name_f);
	imp_calc_f(x,newx,zero_order_info());
	num_f_evals_++;
}

void NLPInterfacePack::NLP::calc_c(const VectorWithOp& x, bool newx) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_.c, "NLP::calc_c()", name_c);
	imp_calc_c(x,newx,zero_order_info());
	num_c_evals_++;
}

void NLPInterfacePack::NLP::calc_h(const VectorWithOp& x, bool newx) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_.h, "NLP::calc_h()", name_h);
	imp_calc_h(x,newx,zero_order_info());
	num_h_evals_++;
}

void NLPInterfacePack::NLP::report_final_solution(
	const VectorWithOp&    x
	,const VectorWithOp*   lambda
	,const VectorWithOp*   lambdaI
	,const VectorWithOp*   nu
	,bool                  optimal
	) const
{
	// The default behavior is just to ignore this!
}

NLPInterfacePack::size_type NLPInterfacePack::NLP::num_f_evals() const
{
	return num_f_evals_;
}

NLPInterfacePack::size_type NLPInterfacePack::NLP::num_c_evals() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return num_c_evals_;
}

NLPInterfacePack::size_type NLPInterfacePack::NLP::num_h_evals() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return num_h_evals_;
}
