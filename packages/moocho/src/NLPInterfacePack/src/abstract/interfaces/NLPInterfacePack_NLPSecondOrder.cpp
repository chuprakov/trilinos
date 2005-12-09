// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrder.cpp
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

#include "NLPInterfacePack_NLPSecondOrder.hpp"
#include "Teuchos_TestForException.hpp"

namespace {
	const char name_HL[] = "HL";
}

namespace NLPInterfacePack {

// constructors

NLPSecondOrder::NLPSecondOrder()
	: HL_(NULL)
{}


void NLPSecondOrder::initialize(bool test_setup) {
	num_HL_evals_ = 0;
	NLPFirstOrder::initialize(test_setup);
}

// <<std aggr>> members for HL

void NLPSecondOrder::set_HL(MatrixSymOp* HL)
{
	HL_ = HL;
}

MatrixSymOp* NLPSecondOrder::get_HL()
{
	return StandardCompositionRelationshipsPack::get_role_name(HL_, false, name_HL);
}

MatrixSymOp& NLPSecondOrder::HL()
{
	return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

const MatrixSymOp& NLPSecondOrder::HL() const
{
	return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

void NLPSecondOrder::unset_quantities()
{
	NLPFirstOrder::unset_quantities();
	HL_ = NULL;
}

// calculations

void NLPSecondOrder::calc_HL(
	const Vector& x, const Vector* lambda, bool newpoint
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( lambda  && this->m()  == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(HL_, "NLP::calc_HL()", name_HL);
	imp_calc_HL(x,lambda,newpoint,second_order_info());
}

size_type NLPSecondOrder::num_HL_evals() const
{
	return num_HL_evals_;
}

} // namespace NLPInterfacePack
