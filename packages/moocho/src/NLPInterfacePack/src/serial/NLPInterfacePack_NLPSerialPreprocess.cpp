// ///////////////////////////////////////////////////////////////////////
// NLPSerialPreprocess.cpp
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

#include <assert.h>

#include <iostream> // Debug only
#include "LinAlgPack/include/PermOut.h"

#include <algorithm>
#include <sstream>
#include <limits>
#include <stdio.h>
#include <fstream>

#include "NLPInterfacePack/include/NLPSerialPreprocess.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/PermutationSerial.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/IVector.h"
#include "LinAlgPack/include/PermVecMat.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"
#include "AbstractFactoryStd.h"
#include "dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
}

namespace NLPInterfacePack {

// NLPSerialPreprocess

// Static public members

value_type
NLPSerialPreprocess::fixed_var_mult()
{
	return std::numeric_limits<LinAlgPack::Vector::value_type>::max()-100; // Don't know what to use?
}

// Constructors / nitializers

NLPSerialPreprocess::NLPSerialPreprocess(
	bool  convert_inequ_to_equ
	)
	:convert_inequ_to_equ_(convert_inequ_to_equ)
	,initialized_(false)
	,force_xinit_in_bounds_(true)
	,scale_f_(1.0)
{}

void NLPSerialPreprocess::set_convert_inequ_to_equ( bool convert_inequ_to_equ )
{
	convert_inequ_to_equ_ = convert_inequ_to_equ;
	initialized_          = false;
}

// Overridden public members from NLP

void NLPSerialPreprocess::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
	force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool NLPSerialPreprocess::force_xinit_in_bounds() const
{
	return force_xinit_in_bounds_;
}

void NLPSerialPreprocess::initialize(bool test_setup)
{
	namespace mmp = MemMngPack;

	const value_type inf_bnd = NLP::infinite_bound();

	if( initialized_  && !imp_nlp_has_changed() ) {
		// The subclass NLP has not changed so we can just
		// slip this preprocessing.
		NLPObjGradient::initialize(test_setup);
		return;
	}

	// Get the dimensions of the original problem

	n_orig_  = imp_n_orig();
	m_orig_  = imp_m_orig();   // This may be zero!
	mI_orig_ = imp_mI_orig();  // This may be zero!
	
	// Get the dimensions of the full problem

	n_full_  = n_orig_  + ( convert_inequ_to_equ_ ? mI_orig_ : 0         );
	m_full_  = m_orig_  + ( convert_inequ_to_equ_ ? mI_orig_ : 0         );
	mI_full_ =            ( convert_inequ_to_equ_ ? 0        : mI_orig_  );

	// Initialize the storage for the intermediate quanities
	
	xinit_full_.resize(n_full_);
	xl_full_.resize(n_full_);
	xu_full_.resize(n_full_);
	x_full_.resize(n_full_);
	c_orig_.resize(m_orig_);
	h_orig_.resize(mI_orig_);
	Gf_full_.resize(n_full_);
	var_full_to_fixed_.resize(n_full_);
	equ_perm_.resize(m_full_);
	inv_equ_perm_.resize(m_full_);
	space_c_.initialize(m_full_);
	space_h_.initialize(mI_full_);
	factory_P_var_   = mmp::rcp( new mmp::AbstractFactoryStd<Permutation,PermutationSerial>() );
	factory_P_equ_   = mmp::rcp( new mmp::AbstractFactoryStd<Permutation,PermutationSerial>() );
	factory_P_inequ_ = mmp::rcp( new mmp::AbstractFactoryStd<Permutation,PermutationSerial>() );

	// Intialize xinit_full_, xl_full_ and xu_full_ for the initial point which will set the
	// fixed elements which will not change during the optimization.
	xinit_full_(1,n_orig_)  = imp_xinit_orig();
	xl_full_(1,n_orig_)     = imp_xl_orig();
	xu_full_(1,n_orig_)     = imp_xu_orig();
	if( n_full_ > n_orig_ ) { // Include slack varaibles
		xinit_full_(n_orig_+1,n_full_)  = 0.0;
		xl_full_(n_orig_+1,n_full_)     = imp_hl_orig();
		xu_full_(n_orig_+1,n_full_)     = imp_hu_orig();
	}

	const bool has_var_bounds = imp_has_var_bounds() || n_full_ > n_orig_;

	// Force the initial point in bounds if it is not.
	if( force_xinit_in_bounds() && has_var_bounds ) {
		AbstractLinAlgPack::force_in_bounds(
			VectorWithOpMutableDense( xl_full_(), mmp::null )
			,VectorWithOpMutableDense( xu_full_(), mmp::null )
			,&VectorWithOpMutableDense( x_full_(), mmp::null )
			);
	}
	
	// Determine which variables are fixed by bounds!
	size_type
		xl_nz     = 0,
		xu_nz     = 0,
		num_bnd_x = 0;
	if( has_var_bounds ) {
		// Determine which variables are fixed by bounds and
		// adjust the bounds if needed.
		Vector::iterator
			xl_full		= xl_full_.begin(),
			xu_full		= xu_full_.begin();
		n_ = 0;
		size_type num_fixed = 0;
		for(int i = 1; i <= n_full_; ++i, ++xl_full, ++xu_full) {
			THROW_EXCEPTION(
				*xl_full > *xu_full, InconsistantBounds
				,"NLPSerialPreprocess::initialize() : Error, Inconsistant bounds: xl_full("
				<< i << ") > xu_full(" << i << ")" ); 
			if(*xl_full == *xu_full) {
				//
				// Fixed between bounds
				//
				var_full_to_fixed_(n_full_ - num_fixed) = i;
				num_fixed++;
			}
			else {
				//
				// Not Fixed between bounds
				//
				// Adjust the bounds if needed
				*xl_full = *xl_full < -inf_bnd ? -inf_bnd : *xl_full;
				*xu_full = *xu_full > +inf_bnd ? +inf_bnd : *xu_full;
				//
				n_++;
				var_full_to_fixed_(n_) = i;
				// Check if xl is bounded
				if(	*xl_full != -inf_bnd )
					++xl_nz;
				// Check if xu is bounded
				if(	*xu_full != inf_bnd )
					++xu_nz;
				if( *xl_full != -inf_bnd || *xu_full != inf_bnd )
					++num_bnd_x;
			}
		}
	}
	else {
		// None of the variables are fixed by bounds because there are no bounds
		n_ = n_full_;
		LinAlgPack::identity_perm( &var_full_to_fixed_ );
	}

//	std::cerr << "n_ =" << n_ << std::endl;
//	std::cerr << "var_full_to_fixed_ =\n" << var_full_to_fixed_;
	
	num_bounded_x_ = num_bnd_x;

	// Validate that we still have a solvable problem
	THROW_EXCEPTION(
		n_ < m_full_, InvalidInitialization
		,"NLPSerialPreprocess::initialize() : Error, after removing fixed "
		<< "variables, n = " << n_ << " < m = " << m_full_
		<< ", and the NLP is over determined and can not be solved!" );

	// Initialize inverse of var_full_to_fixed_
	LinAlgPack::inv_perm( var_full_to_fixed_, &inv_var_full_to_fixed_ );

//	std::cerr << "inv_var_full_to_fixed_ =\n"  << inv_var_full_to_fixed_;

	var_perm_.resize(n_);
	space_x_.initialize(n_);
	
	// Resize xinit, xl, xu, hl and hu
	xinit_.initialize(n_);
	if(has_var_bounds) {
		xl_.initialize(n_);
		xu_.initialize(n_);
	}
	if(mI_full_) {
		hl_.initialize(mI_full_);
		hu_.initialize(mI_full_);
	}

	if( m_full_ ) {
		// Get the first basis
		if( !nlp_selects_basis() ) {
			// The NLP is not selecting the first basis so set to the initial basis to
			// the indentity permutations and assume full column rank for Gc.
			LinAlgPack::identity_perm(&var_perm_);
//			std::cerr << "var_perm_ =\n" << var_perm_;
			LinAlgPack::identity_perm(&equ_perm_);
//			std::cerr << "equ_perm_ =\n" << equ_perm_;
			LinAlgPack::identity_perm(&inv_equ_perm_);
//			std::cerr << "inv_equ_perm_ =\n" << inv_equ_perm_;
			r_ = m_full_;
			var_from_full( xinit_full_().begin(), xinit_.set_vec().begin() );
			if(has_var_bounds) {
				var_from_full( xl_full_().begin(), xl_.set_vec().begin() );
				var_from_full( xu_full_().begin(), xu_.set_vec().begin() );
				do_force_xinit_in_bounds();
			}
		}
		else {
			// The nlp subclass is selecting the first basis.
			
			// make intialized_ true temporaraly so you can call get_next_basis()
			// and assert_initialized() called in it will not throw an exception.
			initialized_ = true;
			
			try {
				size_type rank;
				const bool 
 					get_next_basis_return = get_next_basis_remove_fixed( &var_perm_, &equ_perm_, &rank );
				THROW_EXCEPTION(
					!get_next_basis_return, std::logic_error
					,"NLPSerialPreprocess::initialize():  "
					" If nlp_selects_basis() is true then imp_get_next_basis() "
					" must return true for the first call" );
				assert_and_set_basis( var_perm_, equ_perm_, rank );
//				std::cerr << "var_perm_ =\n" << var_perm_;
//				std::cerr << "equ_perm_ =\n" << equ_perm_;
			}
			catch(...) {
				// In case an exception was thrown I don't want to leave #this#
				// in an inconsistant state.
				initialized_ = false;
				throw;
			}
			
			initialized_ = false;	// resize to false to continue initialization
		}
	}

//	std::cerr << "n_full_ = " << n_full_ << std::endl;
//	std::cerr << "n_ = " << n_ << std::endl;
//	std::cerr << "var_full_to_fixed_ =\n" << var_full_to_fixed_;
//	std::cerr << "inv_var_full_to_fixed_ =\n"  << inv_var_full_to_fixed_;
//	std::cerr << "var_perm_ =\n" << var_perm_;
//	std::cerr << "equ_perm_ =\n" << equ_perm_;

	// If you get here then the initialization went Ok.
	NLPObjGradient::initialize(test_setup);
	initialized_ = true;
}

bool NLPSerialPreprocess::is_initialized() const
{
	return initialized_;
}

size_type NLPSerialPreprocess::n() const 
{
	assert_initialized();
	return n_;
}

size_type NLPSerialPreprocess::m() const 
{	
	assert_initialized(); 
	return m_full_;
}

size_type NLPSerialPreprocess::mI() const 
{	
	assert_initialized(); 
	return mI_full_;
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_x() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(&space_x_,false);
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_c() const
{
	namespace mmp = MemMngPack;
	return ( m_full_ ? mmp::rcp(&space_c_,false) : mmp::null );
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_h() const
{
	namespace mmp = MemMngPack;
	return ( mI_full_ ? mmp::rcp(&space_h_,false) : mmp::null );
}

size_type NLPSerialPreprocess::num_bounded_x() const 
{
	return num_bounded_x_;
}

const VectorWithOp& NLPSerialPreprocess::xl() const 
{
	assert_initialized();
	THROW_EXCEPTION(
		num_bounded_x_ == 0, NoBounds
		,"NLPSerialPreprocess::xl() : Error, there are no variable bounds!"
		);
	return xl_;
}

const VectorWithOp& NLPSerialPreprocess::xu() const 
{
	assert_initialized();
	THROW_EXCEPTION(
		num_bounded_x_ == 0, NoBounds
		,"NLPSerialPreprocess::xu() : Error, there are no variable bounds!"
		);
	return xu_;
}

const VectorWithOp& NLPSerialPreprocess::hl() const 
{
	assert_initialized();
	return hl_;
}

const VectorWithOp& NLPSerialPreprocess::hu() const 
{
	assert_initialized();
	return hu_;
}

const VectorWithOp& NLPSerialPreprocess::xinit() const 
{
	assert_initialized();
	return xinit_;
}

void NLPSerialPreprocess::get_init_lagrange_mult(
	VectorWithOpMutable*   lambda
	,VectorWithOpMutable*  lambdaI
	,VectorWithOpMutable*  nu
	) const
{
	// ToDo: Get subclass to define what these are!
	if(lambda)
		*lambda   = 0.0;
	if(lambdaI)
		*lambdaI = 0.0;
	if(nu)
		*nu      = 0.0;
}

void NLPSerialPreprocess::scale_f( value_type scale_f )
{
	scale_f_ = scale_f;
}

value_type NLPSerialPreprocess::scale_f() const
{
	return scale_f_;
}

void NLPSerialPreprocess::report_final_solution(
	const VectorWithOp&    x
	,const VectorWithOp*   lambda
	,const VectorWithOp*   lambdaI
	,const VectorWithOp*   nu
	,bool                  is_optimal
	) const
{
	// set x_full
	VectorDenseEncap  x_d(x);
	Vector x_full(n_full_);
	x_full = x_full_;	// set any fixed components (as well as the others at first)
	var_to_full( x_d().begin(), x_full().begin() );	// set the nonfixed components
	// set lambda_full
	Vector lambda_full;
	if( lambda ) {
		VectorDenseEncap lambda_d(*lambda);
		VectorSlice      lambda = lambda_d();
		lambda_full.resize(m_full_);
		for(size_type j = 1; j <= m_full_; j++)
			lambda_full(equ_perm_(j)) = lambda(j);
	}
	// set lambdaI_full
	Vector lambdaI_full;
	if(lambdaI) {
		VectorDenseEncap lambdaI_d(*lambdaI);
		lambdaI_full = lambdaI_d();
	}
	// set nu_full
	Vector nu_full(n_full_);
	if(nu) {
		nu_full = 0.0; // We don't give lagrange multipliers for fixed varaibles!
		// ToDo: Define a special constrant for multiplier values for fixed variables 
		VectorDenseEncap nu_d(*nu);
		var_to_full( nu_d().begin(), nu_full().begin() );	// set the nonfixed components
	}
	// Report the final solution
	VectorSlice
		lambda_orig   = lambda && m_orig_ ? lambda_full(1,m_orig_) : VectorSlice(),
		lambdaI_orig  = ( lambdaI
						  ? lambdaI_full()
						  : ( lambda && m_full_ > m_orig_ 
							  ? lambda_full(m_orig_+1,m_full_)
							  : VectorSlice() ) ),
		nu_orig       = nu ? nu_full(1,n_orig_) : VectorSlice();
	imp_report_orig_final_solution(
		x_full
		,lambda_orig.dim()  ? &lambda_orig  : NULL
		,lambdaI_orig.dim() ? &lambdaI_orig : NULL
		,nu_orig.dim()      ? &nu_orig      : NULL
		,is_optimal
		);
}

// Overridden public members from NLPVarReductPerm

const NLPVarReductPerm::perm_fcty_ptr_t
NLPSerialPreprocess::factory_P_var() const
{
	return factory_P_var_;
}

const NLPVarReductPerm::perm_fcty_ptr_t
NLPSerialPreprocess::factory_P_equ() const
{
	return factory_P_equ_;
}

const NLPVarReductPerm::perm_fcty_ptr_t
NLPSerialPreprocess::factory_P_inequ() const
{
	return factory_P_inequ_;
}

Range1D NLPSerialPreprocess::var_dep() const
{
	assert_initialized();
	return r_ ? Range1D(1,r_) : Range1D::Invalid;
}

Range1D NLPSerialPreprocess::var_indep() const
{
	assert_initialized();
	return Range1D(r_+1,n_);
}

Range1D NLPSerialPreprocess::equ_decomp() const
{
	assert_initialized();
	return r_ ? Range1D(1,r_) : Range1D::Invalid;
}

Range1D NLPSerialPreprocess::equ_undecomp() const
{
	assert_initialized();
	return r_ < m_full_ ? Range1D(r_+1,m_full_) : Range1D::Invalid;
}

Range1D NLPSerialPreprocess::inequ_decomp() const
{
	assert_initialized();
	return Range1D::Invalid; // Decomposed inequalities are not allowed yet!
}

Range1D NLPSerialPreprocess::inequ_undecomp() const
{
	assert_initialized();
	return mI_full_ ? Range1D(1,mI_full_) : Range1D::Invalid;
}

bool NLPSerialPreprocess::get_next_basis(
	Permutation*  P_var,   Range1D* var_dep
	,Permutation* P_equ,   Range1D* equ_decomp
	,Permutation* P_inequ, Range1D* inequ_decomp
	)
{
	assert_initialized();
	size_type rank = 0;
	const bool 
		get_next_basis_return = get_next_basis_remove_fixed( &var_perm_, &equ_perm_, &rank );
	if(get_next_basis_return)
		assert_and_set_basis( var_perm_, equ_perm_, rank );
	else
		return false; // The NLP subclass did not have a new basis to give us!
	this->get_basis(P_var,var_dep,P_equ,equ_decomp,P_inequ,inequ_decomp);
	return true;
}

void NLPSerialPreprocess::set_basis(
	const Permutation   &P_var,   const Range1D  &var_dep
	,const Permutation  *P_equ,   const Range1D  *equ_decomp
	,const Permutation  *P_inequ, const Range1D  *inequ_decomp
	)
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	THROW_EXCEPTION(
		(m_full_ > 0 && (P_equ == NULL || equ_decomp == NULL))
		|| (mI_full_ > 0 && (P_inequ == NULL || inequ_decomp == NULL ))
		,std::invalid_argument
		,"NLPSerialPreprocess::set_basis(...) : Error!" );
	THROW_EXCEPTION(
		mI_full_ > 0 && inequ_decomp->size() > 0
		,InvalidBasis
		,"NLPSerialPreprocess::set_basis(...) : Error, can't handle decomposed inequalities yet!" );
	THROW_EXCEPTION(
		m_full_ > 0 && var_dep.size() != equ_decomp->size()
		,InvalidBasis
		,"NLPSerialPreprocess::set_basis(...) : Error!" );
	// Get the concrete types
	const PermutationSerial
		&P_var_s   = dyn_cast<const PermutationSerial>(P_var),
		*P_equ_s   = m_full_  ? &dyn_cast<const PermutationSerial>(*P_equ)   : NULL,
		*P_inequ_s = mI_full_ ? &dyn_cast<const PermutationSerial>(*P_inequ) : NULL;
	// Get the underlying permutation vectors
	mmp::ref_count_ptr<IVector>
		var_perm   = mmp::rcp_const_cast<IVector>(P_var_s.perm()),
		equ_perm   = ( m_full_
					   ? mmp::rcp_const_cast<IVector>(P_equ_s->perm())
					   : mmp::null ),
		inequ_perm = ( mI_full_
					   ? mmp::rcp_const_cast<IVector>(P_inequ_s->perm())
					   : mmp::null );
	THROW_EXCEPTION(
		(m_full_ > 0 && equ_perm.get() == NULL)
		,std::invalid_argument
		,"NLPSerialPreprocess::set_basis(...) : Error, P_equ is not initialized properly!" );
	THROW_EXCEPTION(
		(mI_full_ > 0 && inequ_perm.get() == NULL)
		,std::invalid_argument
		,"NLPSerialPreprocess::set_basis(...) : Error, P_inequ is not initialized properly!" );
	// Set the basis
	assert_and_set_basis( *var_perm, *equ_perm, var_dep.size() );
}

void NLPSerialPreprocess::get_basis(
	Permutation*  P_var,   Range1D* var_dep
	,Permutation* P_equ,   Range1D* equ_decomp
	,Permutation* P_inequ, Range1D* inequ_decomp
	) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	assert_initialized();
	THROW_EXCEPTION(
		P_var == NULL || var_dep == NULL
		|| (m_full_ > 0 && (P_equ == NULL || equ_decomp == NULL))
		|| (mI_full_ > 0 && (P_inequ == NULL || inequ_decomp == NULL ))
		,std::invalid_argument
		,"NLPSerialPreprocess::get_basis(...) : Error!" );
	// Get the concrete types
	PermutationSerial
		&P_var_s   = dyn_cast<PermutationSerial>(*P_var),
		*P_equ_s   = m_full_  ? &dyn_cast<PermutationSerial>(*P_equ)   : NULL,
		*P_inequ_s = mI_full_ ? &dyn_cast<PermutationSerial>(*P_inequ) : NULL;
	// Get the underlying permutation vectors
	mmp::ref_count_ptr<IVector>
		var_perm   = mmp::rcp_const_cast<IVector>(P_var_s.perm()),
		equ_perm   = ( m_full_
					   ? mmp::rcp_const_cast<IVector>(P_equ_s->perm())
					   : mmp::null ),
		inequ_perm = ( mI_full_
					   ? mmp::rcp_const_cast<IVector>(P_inequ_s->perm())
					   : mmp::null );
	// Allocate permutation vectors if none allocated yet or someone else has reference to them
	if( var_perm.get() == NULL || var_perm.count() > 2 ) // P_var reference and my reference
		var_perm = mmp::rcp( new IVector(n_) );
	if( m_full_ && ( equ_perm.get() == NULL || equ_perm.count() > 2 ) ) // P_equ reference and my reference
		equ_perm = mmp::rcp( new IVector(m_full_) );
	if( mI_full_ && ( inequ_perm.get() == NULL || inequ_perm.count() > 2 ) ) // P_inequ reference and my reference
		inequ_perm = mmp::rcp( new IVector(mI_full_) );
	// Copy the basis selection
	(*var_perm)   = var_perm_;
	(*var_dep)    = Range1D(1,r_);
	if(m_full_) {
		(*equ_perm)   = equ_perm_;
		(*equ_decomp) = Range1D(1,r_);
	}
	if(mI_full_) {
		LinAlgPack::identity_perm( inequ_perm.get() ); // No decomposed inequalities right now!
		*inequ_decomp = Range1D::Invalid;
	}
	// Reinitialize the Permutation arguments.
	P_var_s.initialize( var_perm, mmp::null, true );  // Allocate the inverse permuation as well!
	if(m_full_)
		P_equ_s->initialize( equ_perm, mmp::null, true );
	if(mI_full_)
		P_inequ_s->initialize( inequ_perm, mmp::null, true );
}

// Overridden protected members from NLP

void NLPSerialPreprocess::imp_calc_f(
	const VectorWithOp      &x
	,bool                   newx
	,const ZeroOrderInfo    &zero_order_info
	) const
{
	assert_initialized();
	VectorDenseEncap          x_d(x);
	set_x_full( x_d(), newx, &x_full_() );
	imp_calc_f_orig( x_full(), newx, zero_order_orig_info() );
	*zero_order_info.f = scale_f_ * f_orig_;
}

void NLPSerialPreprocess::imp_calc_c(
	const VectorWithOp      &x
	,bool                   newx
	,const ZeroOrderInfo    &zero_order_info
	) const
{
	assert_initialized();
	VectorDenseEncap          x_d(x);
	set_x_full( x_d(), newx, &x_full_() );
	if( m_orig_ )
		imp_calc_c_orig( x_full(), newx, zero_order_orig_info() );
	if( mI_orig_ && convert_inequ_to_equ_ )
		imp_calc_h_orig( x_full(), newx, zero_order_orig_info() );
	VectorDenseMutableEncap  c_d(*zero_order_info.c);
	equ_from_full(
		m_orig_                            ? c_orig_()                     : VectorSlice()
		,mI_orig_ && convert_inequ_to_equ_ ? h_orig_()                     : VectorSlice()
		,mI_orig_ && convert_inequ_to_equ_ ? x_full()(n_orig_+1,n_full_)   : VectorSlice() // s_orig
		,&c_d()
		);
}

void NLPSerialPreprocess::imp_calc_h(
	const VectorWithOp      &x
	,bool                   newx
	,const ZeroOrderInfo    &zero_order_info
	) const
{
	// If this function gets called then this->mI() > 0 must be true
	// which means that convert_inequ_to_equ must be false!
	assert_initialized();
	VectorDenseEncap          x_d(x);
	set_x_full( x_d(), newx, &x_full_() );
	imp_calc_h_orig( x_full(), newx, zero_order_orig_info() );
	VectorDenseMutableEncap  h_d(*zero_order_info.h);
	h_d() = h_orig_(); // Nothing fancy right now
}

// Overridden protected members from NLPObjGradient

void NLPSerialPreprocess::imp_calc_Gf(
	const VectorWithOp      &x
	,bool                   newx
	,const ObjGradInfo      &obj_grad_info
	) const
{
	using LinAlgPack::Vt_S;
	assert_initialized();
	VectorDenseEncap          x_d(x);
	set_x_full( x_d(), newx, &x_full_() );
	if( n_full_ > n_orig_ ) Gf_full_(n_orig_+1,n_full_) = 0.0; // Initialize terms for slacks to zero!
	imp_calc_Gf_orig( x_full(), newx, obj_grad_orig_info() );
	VectorDenseMutableEncap  Gf_d(*obj_grad_info.Gf);
	var_from_full( Gf_full_.begin(), Gf_d().begin() );
	Vt_S( &Gf_d(), scale_f_ );
}

// protected members

bool NLPSerialPreprocess::imp_get_next_basis(
	IVector      *var_perm_full
	,IVector     *equ_perm_full
	,size_type   *rank_full
	,size_type   *rank
	)
{
	return false; // default is that the subclass does not select the basis
}

void NLPSerialPreprocess::assert_initialized() const
{
	THROW_EXCEPTION(
		!initialized_, UnInitialized
		,"NLPSerialPreprocess : The nlp has not been initialized yet" );
}

void NLPSerialPreprocess::set_x_full(
	const VectorSlice& x, bool newx
	,VectorSlice* x_full
	) const
{
	LinAlgPack::assert_vs_sizes(x.dim(),n_);
	if(newx)
		var_to_full(x.begin(), x_full->begin());
}

void NLPSerialPreprocess::var_from_full(
	VectorSlice::const_iterator   vec_full
	,VectorSlice::iterator        vec
	) const
{
//	std::cout << "\nvar_from_full(...) : ";
	for(size_type i = 1; i <= n_; i++) {
		*vec++ = vec_full[var_full_to_fixed_(var_perm_(i)) - 1];
//		std::cout
//			<< "\ni = " << i
//			<< "\nvar_perm_(i) = " << var_perm_(i)
//			<< "\nvar_full_to_fixed_(var_perm_(i)) = " << var_full_to_fixed_(var_perm_(i))
//			<< "\nvec_full[var_full_to_fixed_(var_perm_(i)) - 1] = " << vec_full[var_full_to_fixed_(var_perm_(i)) - 1]
//			<< "\nvec[i] = " << *(vec-1) << "\n\n";
	}		
}

void NLPSerialPreprocess::var_to_full( VectorSlice::const_iterator vec, VectorSlice::iterator vec_full ) const
{
	for(size_type i = 1; i <= n_; i++)
		vec_full[var_full_to_fixed_(var_perm_(i)) - 1] = *vec++;
}

void NLPSerialPreprocess::equ_from_full(
	const VectorSlice   &c_orig
	,const VectorSlice  &h_orig
	,const VectorSlice  &s_orig
	,VectorSlice        *c_full
	) const
{
	size_type i;
	// c_full = [ c_orig; h_orig - s_orig ]
 	for(i = 1; i <= m_orig_; ++i)
		(*c_full)(inv_equ_perm_(i)) = c_orig(i);
 	for(i = 1; i <= mI_orig_; ++i)
		(*c_full)(inv_equ_perm_(m_orig_+i)) = h_orig(i) - s_orig(i);
}

// private members

bool NLPSerialPreprocess::get_next_basis_remove_fixed(
	IVector* var_perm, IVector* equ_perm, size_type* rank
	)
{
	IVector var_perm_full(n_full_);
	equ_perm->resize(m_full_);
	size_type rank_full, rank_fixed_removed;
	if( imp_get_next_basis( &var_perm_full, equ_perm, &rank_full, &rank_fixed_removed ) ) {
//		std::cerr << "var_perm_full =\n"  << var_perm_full;
//		std::cerr << "equ_perm =\n"  << *equ_perm;
//		std::cerr << "rank_full = "  << rank_full << std::endl;
		//
		// The NLP subclass has another basis to select
		//
		// Translate the basis by removing variables fixed by bounds.
		// This is where it is important that var_perm_full is
		// sorted in assending order for basis and non-basis variables
		//
		// This is a little bit of a tricky algorithm.  We have to
		// essentially loop through the set of basic and non-basic
		// variables, remove fixed variables and adjust the indexes
		// of the remaining variables.  Since the set of indexes in
		// the basic and non-basis sets are sorted, this will not
		// be too bad of an algorithm.

		// Iterator for the fixed variables that we are to remove
		IVector::const_iterator     fixed_itr     = var_full_to_fixed_.begin() + n_;
		IVector::const_iterator     fixed_end     = var_full_to_fixed_.end();

		// Iterator for the basis variables
		IVector::iterator           basic_itr  = var_perm_full.begin();
		IVector::iterator           basic_end  = basic_itr + rank_full;

		// Iterator for the non-basis variables
		IVector::iterator           nonbasic_itr  = basic_end;
		IVector::iterator           nonbasic_end  = var_perm_full.end();
		
		// Count the number of fixed basic and non-basic variables
		index_type
			count_fixed          = 0,
			count_basic_fixed    = 0,
			count_nonbasic_fixed = 0;

		// Loop through all of the fixed variables and remove and compact
		for( ; fixed_itr != fixed_end; ++fixed_itr ) {
			const index_type
				next_fixed = ( fixed_itr +1 != fixed_end ? *(fixed_itr+1) : n_full_+1);
			// Bring the basic and nonbasic variables up to this fixed variable
			for( ; *basic_itr < *fixed_itr; ++basic_itr )
				*(basic_itr - count_basic_fixed) = *basic_itr - count_fixed;
			for( ; *nonbasic_itr < *fixed_itr; ++nonbasic_itr )
				*(nonbasic_itr - count_nonbasic_fixed) = *nonbasic_itr - count_fixed;
			// Update the count of the fixed variables
			if( *basic_itr == *fixed_itr ) {
				++count_basic_fixed;
				++basic_itr;
			}
			else {
				assert(*nonbasic_itr == *fixed_itr); // If basic was not fixed then nonbasic better be!
				++count_nonbasic_fixed;
				++nonbasic_itr;

			}
			++count_fixed;
			// Now update the indexes until the next fixed variable
			for( ; *basic_itr < next_fixed; ++basic_itr )
				*(basic_itr - count_basic_fixed) = *basic_itr - count_fixed;
			for( ; *nonbasic_itr < next_fixed; ++nonbasic_itr )
				*(nonbasic_itr - count_nonbasic_fixed) = *nonbasic_itr - count_fixed;
		}
		assert(count_fixed == n_full_ - n_); // Basic check

		var_perm->resize(n_);
		std::copy(
			var_perm_full.begin()
			,var_perm_full.begin() + rank_fixed_removed
			,var_perm->begin()
			);
		std::copy(
			var_perm_full.begin() + rank_full
			,var_perm_full.begin() + rank_full + (n_-rank_fixed_removed)
			,var_perm->begin() + rank_fixed_removed
			);
		*rank = rank_fixed_removed;
		return true;
	}
	else {
	
	// try to find the file giving the basis...
	char ind[17];
	sprintf(ind, "%ld", basis_selection_num_);
	std::string fname = "basis_";
	fname += ind;
	fname += ".sel";
	basis_selection_num_++;

	std::ifstream basis_file(fname.c_str());
	if (basis_file)
		{
		// try to read the basis file
		std::string tags;

		int n;
		basis_file >> tags;
		THROW_EXCEPTION(tags != "n", std::logic_error, "Incorrect basis file format - \"n\" expected, \"" << tags << "\" found.");
		basis_file >> n;
		THROW_EXCEPTION(n == NAN && n == 0, std::logic_error, "Incorrect basis file format - n > 0 \"" << n << "\" found.");

		int m;
		basis_file >> tags;
		THROW_EXCEPTION(tags != "m", std::logic_error, "Incorrect basis file format - \"m\" expected, \"" << tags << "\" found.");
		basis_file >> m;
		THROW_EXCEPTION(m == NAN || m == 0 || m > n , std::logic_error, "Incorrect basis file format - 0 < m <= n expected, \"" << m << "\" found.");
		
		int r;
		basis_file >> tags;
		THROW_EXCEPTION(tags != "rank", std::logic_error, "Incorrect basis file format - \"rank\" expected, \"" << tags << "\" found.");
		basis_file >> r;
		THROW_EXCEPTION(r == NAN || r == 0 || r > m, std::logic_error, "Incorrect basis file format - 0 < rank <= m expected, \"" << r << "\" found.");		
		if (rank)
			{ *rank = r; }

		// var_permutation
		basis_file >> tags;
		THROW_EXCEPTION(tags != "var_perm", std::logic_error, "Incorrect basis file format -\"var_perm\" expected, \"" << tags << "\" found.");
		var_perm->resize(n);
		for (int i=0; i < n; i++)
			{
			int var_index;
			basis_file >> var_index;
			THROW_EXCEPTION(var_index == NAN && var_index >= n, std::logic_error, "Incorrect basis file format for var_perm: 0 <= indice <= n expected, \"" << n << "\" found.");
			(*var_perm)(i) = var_index;
			}

		// eqn_permutation
		basis_file >> tags;
		THROW_EXCEPTION(tags != "equ_perm", std::logic_error, "Incorrect basis file format -\"equ_perm\" expected, \"" << tags << "\" found.");
		equ_perm->resize(r);
		for (int i=0; i < r; i++)
			{
			int equ_index;
			basis_file >> equ_index;
			THROW_EXCEPTION(equ_index == NAN || equ_index >= m, std::logic_error, "Incorrect basis file format for equ_perm: 0 <= indice <= m expected, \"" << n << "\" found.");
			(*equ_perm)(i) = equ_index;
			}

		return true;
		}
	}

	return false;
}

void NLPSerialPreprocess::assert_and_set_basis(
	const IVector& var_perm, const IVector& equ_perm, size_type rank
	)
{
	// Assert that this is a valid basis and set the internal basis.  Also repivot 'xinit', 
	// 'xl', and 'xu'.

	const value_type inf_bnd = NLPSerialPreprocess::infinite_bound();

	// Assert the preconditions
	THROW_EXCEPTION(
		var_perm.size() != n_ || equ_perm.size() != m_full_, std::length_error
		,"NLPSerialPreprocess::set_basis():  The sizes "
		"of the permutation vectors does not match the size of the NLP" );
	THROW_EXCEPTION(
		rank > m_full_, InvalidBasis
		,"NLPSerialPreprocess::set_basis():  The rank "
		"of the basis can not be greater that the number of constraints" );
		
	// Set the basis
	r_ = rank;
	if( &var_perm_ != &var_perm )
		var_perm_ = var_perm;
	if( &equ_perm_ != &equ_perm )
		equ_perm_ = equ_perm;
	LinAlgPack::inv_perm( equ_perm_, &inv_equ_perm_ );

	var_from_full( xinit_full_().begin(), xinit_.set_vec().begin() );
	if(num_bounded_x_) {
		var_from_full( xl_full_().begin(), xl_.set_vec().begin() );
		var_from_full( xu_full_().begin(), xu_.set_vec().begin() );
		do_force_xinit_in_bounds();
	};
}

void NLPSerialPreprocess::assert_bounds_on_variables() const
{
	THROW_EXCEPTION(
		!(imp_has_var_bounds() || n_full_ > n_orig_), NLP::NoBounds
		,"There are no bounds on the variables for this NLP" );
}

void NLPSerialPreprocess::do_force_xinit_in_bounds()
{
	AbstractLinAlgPack::force_in_bounds( xl_, xu_, &xinit_ );
}

} // end namespace NLPInterfacePack
