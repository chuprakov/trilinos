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

#include <math.h>

#include "ExampleNLPBanded.h"
#include "LinAlgPack/include/PermVecMat.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
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
	,bool         nlp_selects_basis
	)
	:is_initialized_(false)
	,nlp_selects_basis_(nlp_selects_basis)
	,basis_selection_was_given_(false)
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
	Gc_full_nz_ = nD_ * ( 1 + 2*(bw_-1) + nI_ ); // Overestimate, ToDo: compute the exact value!
	Gh_full_nz_ = 2*mI_;
	xinit_full_.resize(nD_ + nI_);
	xl_full_.resize(xinit_full_.dim());
	xu_full_.resize(xinit_full_.dim());
	hl_full_.resize(mI_);
	hu_full_.resize(mI_);
	co_full_.resize(nD_ + mU_);
	xinit_full_ = xo;
	xl_full_    = xl;
	xu_full_    = xu;
	hl_full_    = hl;
	hu_full_    = hu;
	co_full_    = 0.0;
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
	inform_new_point(newx);
	*zero_order_info.f = ( 1.0 / 2.0 ) * LinAlgPack::dot( x_full, x_full );
}

void ExampleNLPBanded::imp_calc_c_full(
	const VectorSlice            &x
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	inform_new_point(newx);
	if(c_full_updated_)
		return; // c(x) is already computed in *zero_order_info.c
	Vector
		&c  = *zero_order_info.c;
	const size_type
		num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
		I_remainder = nD_ % nI_;
	size_type j = 0;
	for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
		const size_type num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
		for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
			++j;
			const size_type
				klu = ( j - bw_     >= 0   ? bw_-1 : j-1   ),
				kuu = ( j + bw_ - 1 <= nD_ ? bw_-1 : nD_-j );
			value_type
				&c_j = (c(j) = 10.0 * x(j));
			{for( size_type k = 1; k <= klu; ++k ) {
				c_j -= (3.0 / k) * x(j-k);
			}}
			{for( size_type k = 1; k <= kuu; ++k ) {
				c_j -= (3.0 / k) * x(j+k);
			}}
			const value_type
				term = x(nD_ + q_i) + 1;
			c_j *= (term * term);
			c_j += co_full_(j);
		}
	}
	c_full_updated_ = true;
}

void ExampleNLPBanded::imp_calc_h_full(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	inform_new_point(newx);
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::imp_calc_Gf_full(
	const VectorSlice            &x_full
	,bool                        newx
	,const ObjGradInfoSerial     &obj_grad_info
	) const
{
	inform_new_point(newx);
	*obj_grad_info.Gf = x_full;
}

bool ExampleNLPBanded::imp_get_next_basis(
	IVector      *var_perm
	,IVector     *equ_perm
	,size_type   *rank
	)
{
	if(basis_selection_was_given_)
		return false; // Already gave this basis selection.
	// Select the first nD variables as basis variables which gives
	// a nice banded matrix for the basis matrix C.
	assert(var_perm);
	assert(equ_perm);
	assert(rank);
	var_perm->resize(nD_+nI_);
	equ_perm->resize(nD_+mU_);
	LinAlgPack::identity_perm( var_perm );
	LinAlgPack::identity_perm( equ_perm );
	*rank = nD_;
	basis_selection_was_given_ = true;
	return true;
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
	const VectorSlice& x, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	inform_new_point(newx);
	// Compute c(x) if not already
	this->imp_calc_c_full( x, newx, zero_order_full_info() );
	Vector
		&c = *first_order_expl_info.c; // This must not be NULL!
	// Get references/pointers to data for Gc to be computed/updated.
	index_type
		&Gc_nz = *first_order_expl_info.Gc_nz;
	value_type
		*Gc_val = &(*first_order_expl_info.Gc_val)[0];
	index_type
		*Gc_ivect = ( first_order_expl_info.Gc_ivect
					  ? &(*first_order_expl_info.Gc_ivect)[0] : NULL ),
		*Gc_jvect = ( first_order_expl_info.Gc_jvect
					  ? &(*first_order_expl_info.Gc_jvect)[0] : NULL );
	assert( (Gc_ivect != NULL) == (Gc_jvect != NULL) );
	// Set nonzeros for Gc (in sorted compressed column format w.r.t., i.e. grouped by constraints)
	const size_type
		num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
		I_remainder = nD_ % nI_;
	Gc_nz = 0;
	size_type j = 0;
	for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
		const value_type
			x_q = x(nD_ + q_i),
			x_q_term = (x_q + 1) * (x_q + 1);
		const size_type num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
		for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
			++j;
			const size_type
				klu = ( j - bw_     >= 0   ? bw_-1 : j-1   ),
				kuu = ( j + bw_ - 1 <= nD_ ? bw_-1 : nD_-j );
			//
			{for( index_type k = klu; k >= 1; --k ) {
				++Gc_nz;
				*Gc_val++ = -3.0 / k * x_q_term;
				if(Gc_ivect) {
					*Gc_ivect++ = j - k;
					*Gc_jvect++ = j;
				}
			}}
			//
			++Gc_nz;
			*Gc_val++ = 10.0 * x_q_term;
			if(Gc_ivect) {
				*Gc_ivect++ = j;
				*Gc_jvect++ = j;
			}
			//
			{for( index_type k = 1; k <= kuu; ++k ) {
				++Gc_nz;
				*Gc_val++ = -3.0 / k * x_q_term;
				if(Gc_ivect) {
					*Gc_ivect++ = j + k;
					*Gc_jvect++ = j;
				}
			}}
			//
			++Gc_nz;
			*Gc_val++ = 2.0 * (c(j) - co_full_(j)) / (x_q + 1);
			if(Gc_ivect) {
				*Gc_ivect++ = nD_ + q_i;
				*Gc_jvect++ = j;
			}
		}
	}
}

void ExampleNLPBanded::imp_calc_Gh_full(
	const VectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	inform_new_point(newx);
	assert(0); //  ToDo: Implemenet!
}

// private

void ExampleNLPBanded::assert_is_initialized() const
{
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::inform_new_point(bool newx) const
{
	if(newx) {
		c_full_updated_ = false;
	}
}

}	// end namespace NLPInterfacePack
