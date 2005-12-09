// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrder.cpp
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

#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "Thyra_Range1D.hpp"
#include "Teuchos_TestForException.hpp"

namespace {
	const char name_Gc[] = "Gc";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPFirstOrder::NLPFirstOrder()
	: Gc_(NULL)
{}

void NLPFirstOrder::initialize(bool test_setup) {
	num_Gc_evals_ = 0;
	NLPObjGrad::initialize(test_setup);
}

// BasisSystem

const NLPFirstOrder::basis_sys_ptr_t
NLPFirstOrder::basis_sys() const
{
	return Teuchos::null;
}

// <<std aggr>> members for Gc

void NLPFirstOrder::set_Gc(MatrixOp* Gc)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	Gc_ = Gc;
}

MatrixOp* NLPFirstOrder::get_Gc()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(Gc_, false, name_Gc);
}

MatrixOp& NLPFirstOrder::Gc()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

const MatrixOp& NLPFirstOrder::Gc() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

void NLPFirstOrder::unset_quantities()
{
	NLPObjGrad::unset_quantities();
	Gc_ = NULL;
}

// calculations

void NLPFirstOrder::calc_Gc(const Vector& x, bool newx) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc);
	imp_calc_Gc(x,newx,first_order_info());
	num_Gc_evals_++;
}

size_type NLPFirstOrder::num_Gc_evals() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return num_Gc_evals_;
}

}	// end namespace NLPInterfacePack 
