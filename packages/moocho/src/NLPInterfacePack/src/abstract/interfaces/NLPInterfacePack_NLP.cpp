// ///////////////////////////////////////////////////////////////////////
// NLP.cpp

#include "NLPInterfacePack/include/NLP.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/VectorClass.h"

namespace {
	const char name_f[] = "f";
	const char name_c[] = "c";
} // end namespace

// constructors

NLPInterfacePack::NLP::NLP()
{}

// destructor

NLPInterfacePack::NLP::~NLP()
{}

void NLPInterfacePack::NLP::initialize() {
	num_f_evals_ = num_c_evals_ = 0;
}

void NLPInterfacePack::NLP::get_init_lagrange_mult( Vector* lambda, SpVector* nu ) const
{
	const size_type n = this->n();
	if(lambda) {
		lambda->resize(n);
		(*lambda) = 0.0;
	}
	if(nu) {
		nu->resize(n,n); // No nonzero elements!
	}
}

// <<std comp>> members for f

void NLPInterfacePack::NLP::set_f(value_type* f)
{
	f_c_.f = f;
}

NLPInterfacePack::value_type* NLPInterfacePack::NLP::get_f()
{
	return StandardCompositionRelationshipsPack::get_role_name(f_c_.f, false, name_f);
}

NLPInterfacePack::value_type& NLPInterfacePack::NLP::f()
{
	return StandardCompositionRelationshipsPack::role_name(f_c_.f, false, name_f);
}

const NLPInterfacePack::value_type& NLPInterfacePack::NLP::f() const
{
	return StandardCompositionRelationshipsPack::role_name(f_c_.f, false, name_f);
}

// <<std comp>> members for c

void NLPInterfacePack::NLP::set_c(Vector* c)
{
	f_c_.c = c;
}

NLPInterfacePack::Vector* NLPInterfacePack::NLP::get_c()
{
	return StandardCompositionRelationshipsPack::get_role_name(f_c_.c, false, name_c);
}

NLPInterfacePack::Vector& NLPInterfacePack::NLP::c()
{
	return StandardCompositionRelationshipsPack::role_name(f_c_.c, false, name_c);
}

const NLPInterfacePack::Vector& NLPInterfacePack::NLP::c() const
{
	return StandardCompositionRelationshipsPack::role_name(f_c_.c, false, name_c);
}

// calculations

void NLPInterfacePack::NLP::calc_f(const VectorSlice& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(f_c_.f, "NLP::calc_f()", name_f);
	imp_calc_f(x,newx,zero_order_info());
	num_f_evals_++;
}

void NLPInterfacePack::NLP::calc_c(const VectorSlice& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(f_c_.c, "NLP::calc_c()", name_c);
	imp_calc_c(x,newx,zero_order_info());
	num_c_evals_++;
}

void NLPInterfacePack::NLP::report_final_solution(
	  const VectorSlice&	x
	, const VectorSlice*	lambda
	, const SpVectorSlice*	nu
	, bool					optimal		) const
{
	// The default behavior is just to ignore this!
}

NLPInterfacePack::size_type NLPInterfacePack::NLP::num_f_evals() const
{
	return num_f_evals_;
}

NLPInterfacePack::size_type NLPInterfacePack::NLP::num_c_evals() const
{
	return num_c_evals_;
}
