// ///////////////////////////////////////////////////////////////////////
// NLP.cpp

#include <iostream>

#include "../include/NLP.h"
#include "LinAlgPack/include/VectorClass.h"

// static members

const char NLPInterfacePack::NLP::name_f_[] = "f";
const char NLPInterfacePack::NLP::name_c_[] = "c";

// constructors

NLPInterfacePack::NLP::NLP() : f_(0), owns_f_(false), c_(0), owns_c_(false)
{}

// destructor

NLPInterfacePack::NLP::~NLP()
{
	StandardCompositionRelationshipsPack::destory_container_obj(f_,owns_f_);
	StandardCompositionRelationshipsPack::destory_container_obj(c_,owns_c_);
}

void NLPInterfacePack::NLP::initialize() {
	num_f_evals_ = num_c_evals_ = 0;
}

// <<std comp>> members for f

void NLPInterfacePack::NLP::set_f(value_type* f, bool owns_f)
{
	StandardCompositionRelationshipsPack::set_role_name(f_, owns_f_, name_f_
		, f, owns_f);
}

NLPInterfacePack::value_type* NLPInterfacePack::NLP::get_f()
{
	return StandardCompositionRelationshipsPack::get_role_name(f_, owns_f_, name_f_);
}

void NLPInterfacePack::NLP::set_owns_f(bool owns_f)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(f_, owns_f_, name_f_
		, owns_f);
}

bool NLPInterfacePack::NLP::owns_f()
{
	return StandardCompositionRelationshipsPack::owns_role_name(f_, owns_f_, name_f_);
}

NLPInterfacePack::value_type& NLPInterfacePack::NLP::f()
{
	return StandardCompositionRelationshipsPack::role_name(f_, owns_f_, name_f_);
}

const NLPInterfacePack::value_type& NLPInterfacePack::NLP::f() const
{
	return StandardCompositionRelationshipsPack::role_name(f_, owns_f_, name_f_);
}

// <<std comp>> members for c

void NLPInterfacePack::NLP::set_c(Vector* c, bool owns_c)
{
	StandardCompositionRelationshipsPack::set_role_name(c_, owns_c_, name_c_
		, c, owns_c);
}

NLPInterfacePack::Vector* NLPInterfacePack::NLP::get_c()
{
	return StandardCompositionRelationshipsPack::get_role_name(c_, owns_c_, name_c_);
}

void NLPInterfacePack::NLP::set_owns_c(bool owns_c)
{
	StandardCompositionRelationshipsPack::set_owns_role_name(c_, owns_c_, name_c_
		, owns_c);
}

bool NLPInterfacePack::NLP::owns_c()
{
	return StandardCompositionRelationshipsPack::owns_role_name(c_, owns_c_, name_c_);
}

NLPInterfacePack::Vector& NLPInterfacePack::NLP::c()
{
	return StandardCompositionRelationshipsPack::role_name(c_, owns_c_, name_c_);
}

const NLPInterfacePack::Vector& NLPInterfacePack::NLP::c() const
{
	return StandardCompositionRelationshipsPack::role_name(c_, owns_c_, name_c_);
}

// calculations

void NLPInterfacePack::NLP::calc_f(const Vector& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(f_, "NLP::calc_f()", name_f_);
//	std::cout << "calc_f(...): x(1) = " << x(1) << "\n";
	imp_calc_f(x,newx);
	num_f_evals_++;
}

void NLPInterfacePack::NLP::calc_c(const Vector& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(c_, "NLP::calc_c()", name_c_);
//	std::cout << "calc_c(...): x(1) = " << x(1) << "\n";
	imp_calc_c(x,newx);
	num_c_evals_++;
}

void NLPInterfacePack::NLP::report_final_x(
	  const VectorSlice&	x
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
