// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrderInfo.cpp
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

#include "NLPInterfacePack/src/NLPFirstOrderInfo.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOp.hpp"
#include "Range1D.hpp"
#include "ThrowException.hpp"

namespace {
	const char name_Gc[] = "Gc";
	const char name_Gh[] = "Gh";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPFirstOrderInfo::NLPFirstOrderInfo()
	: Gc_(NULL), Gh_(NULL)
{}

void NLPFirstOrderInfo::initialize(bool test_setup) {
	num_Gc_evals_ = 0;
	num_Gh_evals_ = 0;
	NLPObjGradient::initialize(test_setup);
}

// BasisSystem

const NLPFirstOrderInfo::basis_sys_ptr_t
NLPFirstOrderInfo::basis_sys() const
{
	return MemMngPack::null;
}

// <<std aggr>> members for Gc

void NLPFirstOrderInfo::set_Gc(MatrixWithOp* Gc)
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	Gc_ = Gc;
}

MatrixWithOp* NLPFirstOrderInfo::get_Gc()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(Gc_, false, name_Gc);
}

MatrixWithOp& NLPFirstOrderInfo::Gc()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

const MatrixWithOp& NLPFirstOrderInfo::Gc() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

// <<std aggr>> members for Gh

void NLPFirstOrderInfo::set_Gh(MatrixWithOp* Gh)
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	Gh_ = Gh;
}

MatrixWithOp* NLPFirstOrderInfo::get_Gh()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(Gh_, false, name_Gh);
}

MatrixWithOp& NLPFirstOrderInfo::Gh()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gh_, false, name_Gh);
}

const MatrixWithOp& NLPFirstOrderInfo::Gh() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gh_, false, name_Gh);
}

// calculations

void NLPFirstOrderInfo::calc_Gc(const VectorWithOp& x, bool newx) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc);
	imp_calc_Gc(x,newx,first_order_info());
	num_Gc_evals_++;
}

void NLPFirstOrderInfo::calc_Gh(const VectorWithOp& x, bool newx) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(Gh_, "NLP::calc_Gh()", name_Gh);
	imp_calc_Gh(x,newx,first_order_info());
	num_Gh_evals_++;
}

size_type NLPFirstOrderInfo::num_Gc_evals() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return num_Gc_evals_;
}

size_type NLPFirstOrderInfo::num_Gh_evals() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return num_Gh_evals_;
}

}	// end namespace NLPInterfacePack 
