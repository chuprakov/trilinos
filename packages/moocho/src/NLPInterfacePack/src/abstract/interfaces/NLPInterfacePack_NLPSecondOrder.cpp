// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrderInfo.cpp
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

#include "NLPInterfacePack/include/NLPSecondOrderInfo.h"
#include "ThrowException.h"

namespace {
	const char name_HL[] = "HL";
}

// constructors

NLPInterfacePack::NLPSecondOrderInfo::NLPSecondOrderInfo()
	: HL_(0)
{}


void NLPInterfacePack::NLPSecondOrderInfo::initialize() {
	num_HL_evals_ = 0;
	NLPFirstOrderInfo::initialize();
}

// <<std aggr>> members for HL

void NLPInterfacePack::NLPSecondOrderInfo::set_HL(MatrixSymWithOp* HL)
{
	HL_ = HL;
}

NLPInterfacePack::MatrixSymWithOp* NLPInterfacePack::NLPSecondOrderInfo::get_HL()
{
	return StandardCompositionRelationshipsPack::get_role_name(HL_, false, name_HL);
}

NLPInterfacePack::MatrixSymWithOp& NLPInterfacePack::NLPSecondOrderInfo::HL()
{
	return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

const NLPInterfacePack::MatrixSymWithOp& NLPInterfacePack::NLPSecondOrderInfo::HL() const
{
	return StandardCompositionRelationshipsPack::role_name(HL_, false, name_HL);
}

// calculations

void NLPInterfacePack::NLPSecondOrderInfo::calc_HL(
	const VectorWithOp& x, const VectorWithOp* lambda, const VectorWithOp* lambdaI, bool newpoint
	) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( lambda  && this->m()  == 0, std::logic_error, "" );
	THROW_EXCEPTION( lambdaI && this->mI() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(HL_, "NLP::calc_HL()", name_HL);
	imp_calc_HL(x,lambda,lambdaI,newpoint,second_order_info());
}

NLPInterfacePack::size_type NLPInterfacePack::NLPSecondOrderInfo::num_HL_evals() const
{
	return num_HL_evals_;
}
