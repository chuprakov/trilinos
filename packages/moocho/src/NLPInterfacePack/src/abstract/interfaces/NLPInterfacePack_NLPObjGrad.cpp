// ///////////////////////////////////////////////////////////////////////
// NLPObjGrad.cpp
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

#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"

namespace {
	const char name_Gf[] = "Gf";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPObjGrad::NLPObjGrad()
	: Gf_(NULL)
{}

void NLPObjGrad::initialize(bool test_setup) {
	num_Gf_evals_ = 0;
	NLP::initialize(test_setup);
}

// <<std aggr>> members for Gf

void NLPObjGrad::set_Gf(VectorMutable* Gf)
{
	Gf_ = Gf;
}

AbstractLinAlgPack::VectorMutable* NLPObjGrad::get_Gf()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gf_, false, name_Gf);
}

AbstractLinAlgPack::VectorMutable& NLPObjGrad::Gf()
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

const AbstractLinAlgPack::Vector& NLPObjGrad::Gf() const
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

void NLPObjGrad::unset_quantities()
{
	NLP::unset_quantities();
	Gf_ = NULL;
}

// calculations

void NLPObjGrad::calc_Gf(const Vector& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gf_, "NLP::calc_Gf()", name_Gf);
	imp_calc_Gf(x,newx,obj_grad_info());
	num_Gf_evals_++;
}

size_type NLPObjGrad::num_Gf_evals() const
{
	return num_Gf_evals_;
}

}	// end namespace NLPInterfacePack 
