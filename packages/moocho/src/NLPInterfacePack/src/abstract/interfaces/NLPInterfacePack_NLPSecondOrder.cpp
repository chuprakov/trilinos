// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrderInfo.cpp

#include "../include/NLPSecondOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"

// static memebers


const char NLPInterfacePack::NLPSecondOrderInfo::name_HL_[] = "HL";

// constructors


NLPInterfacePack::NLPSecondOrderInfo::NLPSecondOrderInfo()
	: HL_(0), owns_HL_(false)
{}

// destructor


NLPInterfacePack::NLPSecondOrderInfo::~NLPSecondOrderInfo()
{
	StandardCompositionRelationshipsPack::destory_container_obj(HL_,owns_HL_);
}

void NLPInterfacePack::NLPSecondOrderInfo::initialize() {
	NLPFirstOrderInfo::initialize();
}

// <<std comp>> members for HL


void NLPInterfacePack::NLPSecondOrderInfo::set_HL(MatrixWithOp* HL, bool owns_HL)
{
	StandardCompositionRelationshipsPack::set_role_name(HL_, owns_HL_, name_HL_
		, HL, owns_HL);
}


NLPInterfacePack::MatrixWithOp* NLPInterfacePack::NLPSecondOrderInfo::get_HL()
{
	return StandardCompositionRelationshipsPack::get_role_name(HL_, owns_HL_, name_HL_);
}


void NLPInterfacePack::NLPSecondOrderInfo::set_owns_HL(bool owns_HL)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(HL_, owns_HL_, name_HL_
		, owns_HL);
}


bool NLPInterfacePack::NLPSecondOrderInfo::owns_HL()
{
	return StandardCompositionRelationshipsPack::owns_role_name(HL_, owns_HL_, name_HL_);
}


NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPSecondOrderInfo::HL()
{
	return StandardCompositionRelationshipsPack::role_name(HL_, owns_HL_, name_HL_);
}


const NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPSecondOrderInfo::HL() const
{
	return StandardCompositionRelationshipsPack::role_name(HL_, owns_HL_, name_HL_);
}

// calculations


void NLPInterfacePack::NLPSecondOrderInfo::calc_HL(
	const VectorSlice& x, const VectorSlice& lambda, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(HL_, "NLP::calc_HL()", name_HL_);
	imp_calc_HL(x,lambda,newx);
}
