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

#include "../include/NLPFirstOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorClass.h"

namespace {
	const char name_Gc[] = "Gc";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPFirstOrderInfo::NLPFirstOrderInfo()
	: Gc_(0)
{}

void NLPFirstOrderInfo::initialize() {
	num_Gc_evals_ = 0;
	NLPObjGradient::initialize();
}


// <<std aggr>> members for Gc

void NLPFirstOrderInfo::set_Gc(MatrixWithOp* Gc)
{
	Gc_ = Gc;
}

MatrixWithOp* NLPFirstOrderInfo::get_Gc()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gc_, false, name_Gc);
}

MatrixWithOp& NLPFirstOrderInfo::Gc()
{
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

const MatrixWithOp& NLPFirstOrderInfo::Gc() const
{
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

// calculations

void NLPFirstOrderInfo::calc_Gc(const VectorSlice& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc);
	imp_calc_Gc(x,newx,first_order_info());
	num_Gc_evals_++;
}

size_type NLPFirstOrderInfo::num_Gc_evals() const
{
	return num_Gc_evals_;
}

}	// end namespace NLPInterfacePack 
