// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrderInfo.cpp

#include "../include/NLPSecondOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"

// static memebers


const char NLPInterfacePack::NLPSecondOrderInfo::name_Hf_[] = "Hf";


const char NLPInterfacePack::NLPSecondOrderInfo::name_Hcj_[] = "Hcj";

// constructors


NLPInterfacePack::NLPSecondOrderInfo::NLPSecondOrderInfo()
	: Hf_(0), owns_Hf_(false), Hcj_(0), owns_Hcj_(false)
{}

// destructor


NLPInterfacePack::NLPSecondOrderInfo::~NLPSecondOrderInfo()
{
	StandardCompositionRelationshipsPack::destory_container_obj(Hf_,owns_Hf_);
	StandardCompositionRelationshipsPack::destory_container_obj(Hcj_,owns_Hcj_);
}

void NLPInterfacePack::NLPSecondOrderInfo::initialize() {
	NLPFirstOrderInfo::initialize();
}

// <<std comp>> members for Hf


void NLPInterfacePack::NLPSecondOrderInfo::set_Hf(MatrixWithOp* Hf, bool owns_Hf
	, bool same_struct)
{
	StandardCompositionRelationshipsPack::set_role_name(Hf_, owns_Hf_, name_Hf_
		, Hf, owns_Hf);
}


NLPInterfacePack::MatrixWithOp* NLPInterfacePack::NLPSecondOrderInfo::get_Hf()
{
	return StandardCompositionRelationshipsPack::get_role_name(Hf_, owns_Hf_, name_Hf_);
}


void NLPInterfacePack::NLPSecondOrderInfo::set_owns_Hf(bool owns_Hf)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(Hf_, owns_Hf_, name_Hf_
		, owns_Hf);
}


bool NLPInterfacePack::NLPSecondOrderInfo::owns_Hf()
{
	return StandardCompositionRelationshipsPack::owns_role_name(Hf_, owns_Hf_, name_Hf_);
}


NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPSecondOrderInfo::Hf()
{
	return StandardCompositionRelationshipsPack::role_name(Hf_, owns_Hf_, name_Hf_);
}


const NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPSecondOrderInfo::Hf() const
{
	return StandardCompositionRelationshipsPack::role_name(Hf_, owns_Hf_, name_Hf_);
}

// <<std comp>> members for Hcj


void NLPInterfacePack::NLPSecondOrderInfo::set_Hcj(MatrixWithOp* Hcj, bool owns_Hcj
	, bool same_struct)
{
	StandardCompositionRelationshipsPack::set_role_name(Hcj_, owns_Hcj_, name_Hf_
		, Hcj, owns_Hcj);
}


NLPInterfacePack::MatrixWithOp* NLPInterfacePack::NLPSecondOrderInfo::get_Hcj()
{
	return StandardCompositionRelationshipsPack::get_role_name(Hcj_, owns_Hcj_, name_Hf_);
}


void NLPInterfacePack::NLPSecondOrderInfo::set_owns_Hcj(bool owns_Hcj)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(Hcj_, owns_Hcj_, name_Hf_
		, owns_Hcj);
}


bool NLPInterfacePack::NLPSecondOrderInfo::owns_Hcj()
{
	return StandardCompositionRelationshipsPack::owns_role_name(Hcj_, owns_Hcj_, name_Hf_);
}


NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPSecondOrderInfo::Hcj()
{
	return StandardCompositionRelationshipsPack::role_name(Hcj_, owns_Hcj_, name_Hf_);
}


const NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPSecondOrderInfo::Hcj() const
{
	return StandardCompositionRelationshipsPack::role_name(Hcj_, owns_Hcj_, name_Hf_);
}

// calculations


void NLPInterfacePack::NLPSecondOrderInfo::calc_Hf(
	const Vector& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Hf_, "NLP::calc_Hf()", name_Hf_);
	imp_calc_Hf(x,newx);
}


void NLPInterfacePack::NLPSecondOrderInfo::calc_Hcj(
	const Vector& x, size_type j, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Hcj_, "NLP::calc_Hcj()", name_Hcj_);
	imp_calc_Hcj(x,j,newx);
}