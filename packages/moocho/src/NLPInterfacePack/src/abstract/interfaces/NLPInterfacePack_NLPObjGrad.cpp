// ///////////////////////////////////////////////////////////////////////
// NLPObjGradient.cpp

#include "../include/NLPObjGradient.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorClass.h"

namespace {
	const char name_Gf[] = "Gf";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPObjGradient::NLPObjGradient()
	: Gf_(NULL)
{}

void NLPObjGradient::initialize() {
	num_Gf_evals_ = 0;
	NLP::initialize();
}

// <<std aggr>> members for Gf

void NLPObjGradient::set_Gf(Vector* Gf)
{
	Gf_ = Gf;
}

Vector* NLPObjGradient::get_Gf()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gf_, false, name_Gf);
}


Vector& NLPObjGradient::Gf()
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}


const Vector& NLPObjGradient::Gf() const
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, false, name_Gf);
}

// calculations

void NLPObjGradient::calc_Gf(const VectorSlice& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gf_, "NLP::calc_Gf()", name_Gf);
	imp_calc_Gf(x,newx,obj_grad_info());
	num_Gf_evals_++;
}

size_type NLPObjGradient::num_Gf_evals() const
{
	return num_Gf_evals_;
}

}	// end namespace NLPInterfacePack 
