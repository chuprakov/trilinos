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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include <algorithm>
#include <sstream>
#include <limits>

#include "NLPInterfacePack/include/NLPSerialPreprocess.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/PermutationSerial.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
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

NLPSerialPreprocess::NLPSerialPreprocess()
	:initialized_(false)
	,force_xinit_in_bounds_(true)
	,scale_f_(1.0)
{}
	
// Overridden public members from NLP

void NLPSerialPreprocess::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
	force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool NLPSerialPreprocess::force_xinit_in_bounds() const
{
	return force_xinit_in_bounds_;
}

void NLPSerialPreprocess::initialize()
{
	namespace rcp = MemMngPack;
	namespace afp = MemMngPack;

	const value_type inf_bnd = NLP::infinite_bound();

	if( initialized_  && !imp_nlp_has_changed() ) {
		// The subclass NLP has not changed so we can just
		// slip this preprocessing.
		NLPObjGradient::initialize();
		return;
	}

	// Get the dimension of the full problem

	n_full_  = imp_n_full();
	m_full_  = imp_m_full();   // This may be zero!
	mI_full_ = imp_mI_full();  // This may be zero!
	
	// Initialize the storage for the intermediate quanities
	
	x_full_.resize(n_full_);
	c_full_.resize(m_full_);
	h_full_.resize(mI_full_);
	Gf_full_.resize(n_full_);
	var_full_to_fixed_.resize(n_full_);
	equ_perm_.resize(m_full_);
	space_c_.initialize(m_full_);
	space_h_.initialize(mI_full_);
	factory_P_var_   = rcp::rcp( new afp::AbstractFactoryStd<Permutation,PermutationSerial>() );
	factory_P_equ_   = rcp::rcp( new afp::AbstractFactoryStd<Permutation,PermutationSerial>() );
	factory_P_inequ_ = rcp::rcp( new afp::AbstractFactoryStd<Permutation,PermutationSerial>() );

	// Intialize x_full_ for the initial point which will set the
	// fixed elements which will not change durring the optimization.
	x_full_ = imp_xinit_full();
	// Force the initial point in bounds if it is not.
	if( force_xinit_in_bounds() && imp_has_var_bounds() ) {
		AbstractLinAlgPack::force_in_bounds(
			VectorWithOpMutableDense( imp_xl_full(), rcp::null )
			,VectorWithOpMutableDense( imp_xu_full(), rcp::null )
			,&VectorWithOpMutableDense( x_full_(), rcp::null )
			);
	}

	// Determine which variables are fixed by bounds!
	size_type
		xl_nz     = 0,
		xu_nz     = 0,
		num_bnd_x = 0;
	if( imp_has_var_bounds() ) {
		// Determine which variables are fixed by bounds

		VectorSlice::const_iterator
			xl_full		= imp_xl_full().begin(),
			xu_full		= imp_xu_full().begin();

		n_ = 0;
		size_type num_fixed = 0;
		for(int i = 1; i <= n_full_; ++i, ++xl_full, ++xu_full) {
			THROW_EXCEPTION(
				*xl_full > *xu_full, InconsistantBounds
				,"NLPSerialPreprocess::initialize() : Error, Inconsistant bounds: xl_full("
				<< i << ") > xu_full(" << i << ")" ); 
			if(*xl_full == *xu_full) {
				// Fixed between bounds
				var_full_to_fixed_(n_full_ - num_fixed) = i;
				num_fixed++;
			}
			else {
				// Not Fixed between bounds
				n_++;
				var_full_to_fixed_(n_) = i;
				// Check if xl is bounded
				if(	*xl_full != -inf_bnd ) {
					++xl_nz;
					++num_bnd_x;
				}
				// Check if xu is bounded
				if(	*xu_full != inf_bnd ) {
					++xu_nz;
					++num_bnd_x;
				}
			}
		}
	}
	else {
		// None of the variables are fixed by bounds because there are no bounds
		n_ = n_full_;
		LinAlgPack::identity_perm( &var_full_to_fixed_ );
	}
	
	num_bounded_x_ = num_bnd_x;

	// Validate that we still have a solvable problem
	THROW_EXCEPTION(
		n_ < m_full_, InvalidInitialization
		,"NLPSerialPreprocess::initialize() : Error, after removing fixed "
		<< "variables, n = " << n_ << " < m = " << m_full_
		<< ", and the NLP is over determined and can not be solved!" );

	// Initialize inverse of var_full_to_fixed_
	LinAlgPack::inv_perm( var_full_to_fixed_, &inv_var_full_to_fixed_ );

	var_perm_.resize(n_);
	space_x_.initialize(n_);
	
	// Resize xinit, xl, xu, hl and hu
	xinit_.initialize(n_);
	if(imp_has_var_bounds()) {
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
			LinAlgPack::identity_perm(&equ_perm_);
			r_ = m_full_;
			var_from_full( imp_xinit_full().begin(), xinit_.set_vec().begin() );
			if(imp_has_var_bounds()) {
				var_from_full( imp_xl_full().begin(), xl_.set_vec().begin() );
				var_from_full( imp_xu_full().begin(), xu_.set_vec().begin() );
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

	// If you get here then the initialization went Ok.
	NLPObjGradient::initialize();
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
	namespace rcp = MemMngPack;
	return rcp::rcp(&space_x_,false);
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_c() const
{
	namespace rcp = MemMngPack;
	return ( m_full_ ? rcp::rcp(&space_c_,false) : rcp::null );
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_h() const
{
	namespace rcp = MemMngPack;
	return ( mI_full_ ? rcp::rcp(&space_h_,false) : rcp::null );
}

size_type NLPSerialPreprocess::num_bounded_x() const 
{
	return num_bounded_x_;
}

const VectorWithOp& NLPSerialPreprocess::xl() const 
{
	assert_initialized();
	return xl_;
}

const VectorWithOp& NLPSerialPreprocess::xu() const 
{
	assert_initialized();
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
	SpVector lambdaI_full;
	if(lambdaI) {
		VectorDenseEncap lambdaI_d(*lambdaI);
		Vector lambdaI_full_dense = lambdaI_d();
		compact( lambdaI_full_dense, 0.0, &lambdaI_full );
	}
	// set nu_full
	SpVector nu_full;
	if(nu) {
		VectorDenseEncap nu_d(*nu);
		Vector nu_full_dense(n_full_);
		var_to_full( nu_d().begin(), nu_full_dense().begin() );
		const value_type
			fixed_mult = fixed_var_mult(); 
		for( size_type i = n_+1; i <= n_full_; ++i )
			nu_full_dense(var_full_to_fixed_(i)) = fixed_mult;
		compact( nu_full_dense, 0.0, &nu_full );
	}
	// Report the final solution
	imp_report_full_final_solution(
		x_full()
		,lambda  ? &lambda_full()  : NULL
		,lambdaI ? &lambdaI_full() : NULL
		,nu      ? &nu_full()      : NULL
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
	namespace rcp = MemMngPack;
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
	rcp::ref_count_ptr<IVector>
		var_perm   = rcp::rcp_const_cast<IVector>(P_var_s.perm()),
		equ_perm   = ( m_full_
					   ? rcp::rcp_const_cast<IVector>(P_equ_s->perm())
					   : rcp::null ),
		inequ_perm = ( mI_full_
					   ? rcp::rcp_const_cast<IVector>(P_inequ_s->perm())
					   : rcp::null );
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
	namespace rcp = MemMngPack;
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
	rcp::ref_count_ptr<IVector>
		var_perm   = rcp::rcp_const_cast<IVector>(P_var_s.perm()),
		equ_perm   = ( m_full_
					   ? rcp::rcp_const_cast<IVector>(P_equ_s->perm())
					   : rcp::null ),
		inequ_perm = ( mI_full_
					   ? rcp::rcp_const_cast<IVector>(P_inequ_s->perm())
					   : rcp::null );
	// Allocate permutation vectors if none allocated yet or someone else has reference to them
	if( var_perm.get() == NULL || var_perm.count() > 2 ) // P_var reference and my reference
		var_perm = rcp::rcp( new IVector(n_) );
	if( m_full_ && ( equ_perm.get() == NULL || equ_perm.count() > 2 ) ) // P_equ reference and my reference
		equ_perm = rcp::rcp( new IVector(m_full_) );
	if( mI_full_ && ( inequ_perm.get() == NULL || inequ_perm.count() > 2 ) ) // P_inequ reference and my reference
		inequ_perm = rcp::rcp( new IVector(mI_full_) );
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
	P_var_s.initialize( var_perm, rcp::null, true );  // Allocate the inverse permuation as well!
	if(m_full_)
		P_equ_s->initialize( equ_perm, rcp::null, true );
	if(mI_full_)
		P_inequ_s->initialize( inequ_perm, rcp::null, true );
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
	imp_calc_f_full( x_full(), newx, zero_order_full_info() );
	*zero_order_info.f = scale_f_ * f_full_;
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
	imp_calc_c_full( x_full(), newx, zero_order_full_info() );
	VectorDenseMutableEncap  c_d(*zero_order_info.c);
	equ_from_full( c_full_.begin(), c_d().begin() );	
}

void NLPSerialPreprocess::imp_calc_h(
	const VectorWithOp      &x
	,bool                   newx
	,const ZeroOrderInfo    &zero_order_info
	) const
{
	assert_initialized();
	VectorDenseEncap          x_d(x);
	set_x_full( x_d(), newx, &x_full_() );
	imp_calc_h_full( x_full(), newx, zero_order_full_info() );
	VectorDenseMutableEncap  h_d(*zero_order_info.h);
	h_d() = h_full_(); // Nothing fancy right now
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
	imp_calc_Gf_full( x_full(), newx, obj_grad_full_info() );
	VectorDenseMutableEncap  Gf_d(*obj_grad_info.Gf);
	var_from_full( Gf_full_.begin(), Gf_d().begin() );
	Vt_S( &Gf_d(), scale_f_ );
}

// protected members

bool NLPSerialPreprocess::imp_get_next_basis(
	IVector     *var_perm
	,IVector    *equ_perm
	,size_type  *rank
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
	LinAlgPack::assert_vs_sizes(x.dim(),x_full->dim());
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

void NLPSerialPreprocess::var_to_full(
	VectorSlice::const_iterator    vec
	, VectorSlice::iterator        vec_full
	) const
{
	for(size_type i = 1; i <= n_; i++)
		vec_full[var_full_to_fixed_(var_perm_(i)) - 1] = *vec++;
}

void NLPSerialPreprocess::equ_from_full(
	VectorSlice::const_iterator   vec_full
	,VectorSlice::iterator        vec
	) const
{
	for(size_type i = 1; i <= m_full_; i++)
		*vec++ = vec_full[equ_perm_(i) - 1];
}

// private members

bool NLPSerialPreprocess::get_next_basis_remove_fixed(
	IVector* var_perm, IVector* equ_perm, size_type* rank
	)
{
	IVector var_perm_full(n_full_);
	equ_perm->resize(m_full_);
	size_type rank_full;
	if( imp_get_next_basis( &var_perm_full, equ_perm, &rank_full ) ) {
		//
		// The NLP subclass does have another basis to select
		//
		// Translate the basis by removing variables fixed by bounds.
		// This is where it is important that var_perm is
		// sorted in assending order.

		IVector::const_iterator     itr_fixed     = var_full_to_fixed_.begin() + n_;
		IVector::iterator           itr_var_perm  = var_perm_full.begin();
		size_type                   count_fixed   = 0,
			                        total_fixed   = n_full_ - n_;
		while( count_fixed != total_fixed ) {
			if( *itr_var_perm == *itr_fixed ) {
				// This is a variable fixed by bounds so remove it.
				count_fixed++;
			}
			// Compact and adjust the indice
			*(itr_var_perm - count_fixed) = *itr_var_perm - count_fixed;
			++itr_var_perm;
		}
		// We have found and removed all of the fixed variables from _var_perm
		
		var_perm->resize(n_);
		std::copy( var_perm_full.begin(), var_perm_full.begin() + n_, var_perm->begin() );
		*rank = rank_full;
		return true;
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
	var_from_full( imp_xinit_full().begin(), xinit_.set_vec().begin() );
	if(imp_has_var_bounds()) {
		var_from_full( imp_xl_full().begin(), xl_.set_vec().begin() );
		var_from_full( imp_xu_full().begin(), xu_.set_vec().begin() );
		do_force_xinit_in_bounds();
	};
}

void NLPSerialPreprocess::assert_bounds_on_variables() const
{
	THROW_EXCEPTION(
		!imp_has_var_bounds(), NLP::NoBounds
		,"There are no bounds on the variables for this NLP" );
}

void NLPSerialPreprocess::do_force_xinit_in_bounds()
{
	AbstractLinAlgPack::force_in_bounds( xl_, xu_, &xinit_ );
}

size_type NLPSerialPreprocess::compact(
	const Vector& dense_v, value_type val_exclude
	, SpVector* sp_v
	) const
{
	if( !sp_v->nz() )
		return 0;

	Vector::const_iterator	itr_dense_v	= dense_v.begin();
	SpVector::iterator		itr_sp_v	= sp_v->begin();

	size_type nz = 0;
	for(size_type k = 1; k <= dense_v.dim(); ++k, ++itr_dense_v) {
		if( *itr_dense_v != val_exclude ) {
			itr_sp_v->initialize(k,*itr_dense_v);
			++itr_sp_v;
			nz++;
		}
	}
	return nz;
}

} // end namespace NLPInterfacePack
