// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrderInfo.cpp

#include "../include/NLPFirstOrderInfo.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/VectorClass.h"

// static members

const char NLPInterfacePack::NLPFirstOrderInfo::name_Gf_[] = "Gf";


const char NLPInterfacePack::NLPFirstOrderInfo::name_Gc_[] = "Gc";

// constructors

NLPInterfacePack::NLPFirstOrderInfo::NLPFirstOrderInfo()
	: Gf_(0), owns_Gf_(false), Gc_(0), owns_Gc_(false)
{}

// destructor

NLPInterfacePack::NLPFirstOrderInfo::~NLPFirstOrderInfo()
{
	StandardCompositionRelationshipsPack::destory_container_obj(Gf_,owns_Gf_);
	StandardCompositionRelationshipsPack::destory_container_obj(Gc_,owns_Gc_);
}

void NLPInterfacePack::NLPFirstOrderInfo::initialize() {
	num_Gf_evals_ = num_Gc_evals_ = 0;
	NLP::initialize();
}

// <<std comp>> members for Gf

void NLPInterfacePack::NLPFirstOrderInfo::set_Gf(Vector* Gf, bool owns_Gf)
{
	StandardCompositionRelationshipsPack::set_role_name(Gf_, owns_Gf_, name_Gf_
		, Gf, owns_Gf);
}

NLPInterfacePack::Vector* NLPInterfacePack::NLPFirstOrderInfo::get_Gf()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gf_, owns_Gf_, name_Gf_);
}


void NLPInterfacePack::NLPFirstOrderInfo::set_owns_Gf(bool owns_Gf)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(Gf_, owns_Gf_, name_Gf_
		, owns_Gf);
}


bool NLPInterfacePack::NLPFirstOrderInfo::owns_Gf()
{
	return StandardCompositionRelationshipsPack::owns_role_name(Gf_, owns_Gf_, name_Gf_);
}


NLPInterfacePack::Vector& NLPInterfacePack::NLPFirstOrderInfo::Gf()
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, owns_Gf_, name_Gf_);
}


const NLPInterfacePack::Vector& NLPInterfacePack::NLPFirstOrderInfo::Gf() const
{
	return StandardCompositionRelationshipsPack::role_name(Gf_, owns_Gf_, name_Gf_);
}

// <<std comp>> members for Gc


void NLPInterfacePack::NLPFirstOrderInfo::set_Gc(MatrixWithOp* Gc, bool owns_Gc
	, bool same_struct)
{
	StandardCompositionRelationshipsPack::set_role_name(Gc_, owns_Gc_, name_Gc_
		, Gc, owns_Gc);
}


NLPInterfacePack::MatrixWithOp* NLPInterfacePack::NLPFirstOrderInfo::get_Gc()
{
	return StandardCompositionRelationshipsPack::get_role_name(Gc_, owns_Gc_, name_Gc_);
}


void NLPInterfacePack::NLPFirstOrderInfo::set_owns_Gc(bool owns_Gc)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(Gc_, owns_Gc_, name_Gc_
		, owns_Gc);
}


bool NLPInterfacePack::NLPFirstOrderInfo::owns_Gc()
{
	return StandardCompositionRelationshipsPack::owns_role_name(Gc_, owns_Gc_, name_Gc_);
}


NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPFirstOrderInfo::Gc()
{
	return StandardCompositionRelationshipsPack::role_name(Gc_, owns_Gc_, name_Gc_);
}


const NLPInterfacePack::MatrixWithOp& NLPInterfacePack::NLPFirstOrderInfo::Gc() const
{
	return StandardCompositionRelationshipsPack::role_name(Gc_, owns_Gc_, name_Gc_);
}

// calculations


void NLPInterfacePack::NLPFirstOrderInfo::calc_Gf(const Vector& x
	, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gf_, "NLP::calc_Gf()", name_Gf_);
	imp_calc_Gf(x,newx);
	num_Gf_evals_++;
}


void NLPInterfacePack::NLPFirstOrderInfo::calc_Gc(const Vector& x
	, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc_);
	imp_calc_Gc(x,newx);
	num_Gc_evals_++;
}

NLPInterfacePack::size_type NLPInterfacePack::NLPFirstOrderInfo::num_Gf_evals() const
{
	return num_Gf_evals_;
}

NLPInterfacePack::size_type NLPInterfacePack::NLPFirstOrderInfo::num_Gc_evals() const
{
	return num_Gc_evals_;
}