// //////////////////////////////////////////
// ExampleNLPBanded.cpp
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

#include "ExampleNLPBanded.h"
#include "LinAlgPack/include/VectorOp.h"
#include "ThrowException.h"

namespace NLPInterfacePack {

// Constructors / initializers

ExampleNLPBanded::ExampleNLPBanded(
	size_type     nD
	,size_type    nI
	,size_type    bw
	,size_type    mU
	,size_type    mI
	,value_type   xo
	,value_type   xl
	,value_type   xu
	,value_type   hl
	,value_type   hu
	)
	:is_initialized_(false)
	,nlp_selects_basis_(false)
	,multi_calc_(false)
	,nD_(nD)
	,nI_(nI)
	,bw_(bw)
	,mU_(mU)
	,mI_(mI)
{
#ifdef _DEBUG	
	const char msg_err_head[] = "ExampleNLPBanded::ExampleNLPBanded(...) : Error";
	THROW_EXCEPTION(
		nI > nD, std::invalid_argument
		,msg_err_head<<"!" );
	THROW_EXCEPTION(
		bw < 1 || nD < bw, std::invalid_argument
		,msg_err_head<<"!" );
	THROW_EXCEPTION(
		mU != 0, std::invalid_argument
		,msg_err_head<<", can't handle undecomposed equalities yet!" );
	THROW_EXCEPTION(
		mI != 0, std::invalid_argument
		,msg_err_head<<", can't handle general equalities yet!" );
#endif
	Gc_full_nz_ = bw_ * nD_; // Overestimate, ToDo: compute the exact value!
	Gh_full_nz_ = 2*mI_;
	xinit_full_.resize(nD_ + nI_);
	xl_full_.resize(xinit_full_.dim());
	xu_full_.resize(xinit_full_.dim());
	hl_full_.resize(mI_);
	hu_full_.resize(mI_);
	xinit_full_ = xo;
	xl_full_    = xl;
	xu_full_    = xu;
	hl_full_    = hl;
	hu_full_    = hu;
}

// Overridden public members from NLP

void ExampleNLPBanded::initialize()
{
	if(is_initialized_) {
		NLPSerialPreprocessExplJac::initialize();
		return;
	}
	// Nothing to initialize?
	NLPSerialPreprocessExplJac::initialize();
	is_initialized_ = true;
}

bool ExampleNLPBanded::is_initialized() const
{
	return is_initialized_;
}

value_type ExampleNLPBanded::max_var_bounds_viol() const
{
	return +1e+20; // Functions defined everywhere!
}

void ExampleNLPBanded::set_multi_calc(bool multi_calc) const
{
	multi_calc_ = multi_calc;
}

bool ExampleNLPBanded::multi_calc() const
{
	return multi_calc_;
}

// Overridden from NLPVarReductPerm

bool ExampleNLPBanded::nlp_selects_basis() const
{
	return nlp_selects_basis_;
}

// Overridden protected methods from NLPSerialPreprocess */

bool ExampleNLPBanded::imp_nlp_has_changed() const
{
	return !is_initialized_;
}

size_type ExampleNLPBanded::imp_n_full() const
{
	return nD_ + nI_;
}

size_type ExampleNLPBanded::imp_m_full() const
{
	return nD_ + mU_;
}

size_type ExampleNLPBanded::imp_mI_full() const
{
	return mI_;
}

const VectorSlice ExampleNLPBanded::imp_xinit_full() const
{
	return xinit_full_();
}

bool ExampleNLPBanded::imp_has_var_bounds() const
{
	return false; // ToDo: Add bounds!
}

const VectorSlice ExampleNLPBanded::imp_xl_full() const
{
	return xl_full_();
}

const VectorSlice ExampleNLPBanded::imp_xu_full() const
{
	return xu_full_();
}

const VectorSlice ExampleNLPBanded::imp_hl_full() const
{
	return hl_full_();
}

const VectorSlice ExampleNLPBanded::imp_hu_full() const
{
	return hu_full_();
}

void ExampleNLPBanded::imp_calc_f_full(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	*zero_order_info.f = ( 1.0 / 2.0 ) * LinAlgPack::dot( x_full, x_full );
}

void ExampleNLPBanded::imp_calc_c_full(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::imp_calc_h_full(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::imp_calc_Gf_full(
	const VectorSlice            &x_full
	,bool                        newx
	,const ObjGradInfoSerial     &obj_grad_info
	) const
{
	assert(0); //  ToDo: Implemenet!
}

bool ExampleNLPBanded::imp_get_next_basis(
	IVector      *var_perm
	,IVector     *equ_perm
	,size_type   *rank
	)
{
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::imp_report_full_final_solution(
	const VectorSlice      &x_full
	,const VectorSlice     *lambda_full
	,const SpVectorSlice   *lambdaI_full
	,const SpVectorSlice   *nu_full
	,bool                  optimal
	) const
{
	// ToDo: Do something with the final soltuion?
}

// Overridden protected methods from NLPSerialPreprocessExplJac */

size_type ExampleNLPBanded::imp_Gc_nz_full() const
{
	return Gc_full_nz_;
}

size_type ExampleNLPBanded::imp_Gh_nz_full() const
{
	return Gh_full_nz_;
}

void ExampleNLPBanded::imp_calc_Gc_full(
	const VectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::imp_calc_Gh_full(
	const VectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	assert(0); //  ToDo: Implemenet!
}

// private

void ExampleNLPBanded::assert_is_initialized() const
{
	assert(0); //  ToDo: Implemenet!
}

}	// end namespace NLPInterfacePack
