// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrderInfo.cpp

#include "../include/NLPSecondOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"

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

void NLPInterfacePack::NLPSecondOrderInfo::calc_HL( const VectorSlice& x
	, const VectorSlice& lambda, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(HL_, "NLP::calc_HL()", name_HL);
	imp_calc_HL(x,lambda,newx,second_order_info());
}

NLPInterfacePack::size_type NLPInterfacePack::NLPSecondOrderInfo::num_HL_evals() const
{
	return num_HL_evals_;
}
