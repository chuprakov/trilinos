// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrderInfo.cpp

#include "../include/NLPFirstOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorClass.h"


namespace NLPInterfacePack {

// static members

const char NLPFirstOrderInfo::name_Gf_[] = "Gf";


const char NLPFirstOrderInfo::name_Gc_[] = "Gc";

// constructors

NLPFirstOrderInfo::NLPFirstOrderInfo()
	: Gf_(0), owns_Gf_(false), Gc_(0), owns_Gc_(false)
{}

// destructor

NLPFirstOrderInfo::~NLPFirstOrderInfo()
{
	StandardCompositionRelationshipsPack::destory_container_obj(Gf_,owns_Gf_);
	StandardCompositionRelationshipsPack::destory_container_obj(Gc_,owns_Gc_);
}

void NLPFirstOrderInfo::initialize() {
	num_Gf_evals_ = num_Gc_evals_ = 0;
	NLP::initialize();
}

// <<std comp>> members for Gf

void NLPFirstOrderInfo::set_Gf(Vector* Gf, bool owns_Gf)
{
	StandardCompositionRelationshipsPack::set_role_name(Gf_, owns_Gf_, name_Gf_
		, Gf, owns_Gf);
}

Vector* NLPFirstOrderInfo::get_Gf()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gf_, owns_Gf_, name_Gf_);
}


void NLPFirstOrderInfo::set_owns_Gf(bool owns_Gf)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(Gf_, owns_Gf_, name_Gf_
		, owns_Gf);
}


bool NLPFirstOrderInfo::owns_Gf()
{
	return StandardCompositionRelationshipsPack::owns_role_name(Gf_, owns_Gf_, name_Gf_);
}


Vector& NLPFirstOrderInfo::Gf()
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, owns_Gf_, name_Gf_);
}


const Vector& NLPFirstOrderInfo::Gf() const
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, owns_Gf_, name_Gf_);
}

// <<std comp>> members for Gc


void NLPFirstOrderInfo::set_Gc(MatrixWithOp* Gc, bool owns_Gc
	, bool same_struct)
{
	StandardCompositionRelationshipsPack::set_role_name(Gc_, owns_Gc_, name_Gc_
		, Gc, owns_Gc);
}


MatrixWithOp* NLPFirstOrderInfo::get_Gc()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gc_, owns_Gc_, name_Gc_);
}


void NLPFirstOrderInfo::set_owns_Gc(bool owns_Gc)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(Gc_, owns_Gc_, name_Gc_
		, owns_Gc);
}


bool NLPFirstOrderInfo::owns_Gc()
{
	return StandardCompositionRelationshipsPack::owns_role_name(Gc_, owns_Gc_, name_Gc_);
}


MatrixWithOp& NLPFirstOrderInfo::Gc()
{
	return StandardCompositionRelationshipsPack::role_name(Gc_, owns_Gc_, name_Gc_);
}


const MatrixWithOp& NLPFirstOrderInfo::Gc() const
{
	return StandardCompositionRelationshipsPack::role_name(Gc_, owns_Gc_, name_Gc_);
}

// calculations


void NLPFirstOrderInfo::calc_Gf(const VectorSlice& x
	, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gf_, "NLP::calc_Gf()", name_Gf_);
	imp_calc_Gf(x,newx);
	num_Gf_evals_++;
}

void NLPFirstOrderInfo::calc_Gc(const VectorSlice& x
	, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc_);
	imp_calc_Gc(x,newx);
	num_Gc_evals_++;
}

size_type NLPFirstOrderInfo::num_Gf_evals() const
{
	return num_Gf_evals_;
}

size_type NLPFirstOrderInfo::num_Gc_evals() const
{
	return num_Gc_evals_;
}

}	// end namespace NLPInterfacePack 