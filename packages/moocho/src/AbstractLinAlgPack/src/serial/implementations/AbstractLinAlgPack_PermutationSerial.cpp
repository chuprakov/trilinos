// /////////////////////////////////////////////////////////////////////////////
// PermutationSerial.cpp
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

#include "SparseLinAlgPack/include/PermutationSerial.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "LinAlgPack/include/IVector.h"
#include "LinAlgPack/include/PermVecMat.h"
#include "LinAlgPack/include/PermOut.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

// Constructors / initializers

PermutationSerial::PermutationSerial( size_type dim )
	:space_(dim)
{}

PermutationSerial::PermutationSerial(
	const i_vector_ptr_t      &perm
	,const i_vector_ptr_t     &inv_perm
	,bool                     allocate_missing_perm
	,bool                     check_inv_perm
	)
{
	this->initialize(perm,inv_perm,allocate_missing_perm,check_inv_perm);
}
	
void PermutationSerial::initialize_identity( size_type dim )
{
	namespace rcp = MemMngPack;
	space_.initialize(dim);
	perm_     = rcp::null;
	inv_perm_ = rcp::null;
}

void PermutationSerial::initialize(
	const i_vector_ptr_t      &perm
	,const i_vector_ptr_t     &inv_perm
	,bool                     allocate_missing_perm
	,bool                     check_inv_perm
	)
{
	THROW_EXCEPTION(
		perm.get() == NULL && inv_perm.get() == NULL, std::invalid_argument
		,"PermutationSerial::initialize(...) : Error!" );
	if( perm.get() != NULL && inv_perm.get() != NULL ) {
		THROW_EXCEPTION(
			perm->size() != inv_perm->size(), std::invalid_argument
			,"PermutationSerial::initialize(...) : Error!" );
		if(check_inv_perm) {
			// ToDo: Permform this check!
		}
	}
	space_.initialize( perm.get() ? perm->size() : inv_perm->size() );
	perm_     = perm;
	inv_perm_ = inv_perm;
	if( allocate_missing_perm && perm_.get() == NULL ) {
		MemMngPack::ref_count_ptr<IVector>
			_perm(new IVector(inv_perm_->size()));
		LinAlgPack::inv_perm( *inv_perm_, _perm.get() );
		perm_ = _perm;
	}
	if( allocate_missing_perm && inv_perm_.get() == NULL ) {
		MemMngPack::ref_count_ptr<IVector>
			_inv_perm(new IVector(perm_->size()));
		LinAlgPack::inv_perm( *perm_, _inv_perm.get() );
		inv_perm_ = _inv_perm;
	}
}

// Overridden from Permutation

const VectorSpace& PermutationSerial::space() const
{
	return space_;
}

size_type PermutationSerial::dim() const
{
	return space_.dim();
}

bool PermutationSerial::is_identity() const
{
	return perm_.get() == NULL && inv_perm_.get() == NULL;
}

std::ostream& PermutationSerial::output(std::ostream& out) const
{
	const size_type dim = this->dim();
	out << "Serial " << dim << " x " << dim << " permtutation matrix:\n";
	out << "perm =";
	if( perm_.get() )
		out << "\n" << *perm_;
	else
		out << " NULL\n";
	out << "inv_perm =";
	if( inv_perm_.get() )
		out << "\n" << *inv_perm_;
	else
		out << " NULL\n";
	return out;
}

void PermutationSerial::permute( 
	BLAS_Cpp::Transp          P_trans
	,const VectorWithOp       &x
	,VectorWithOpMutable      *y
	) const
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		y == NULL, std::invalid_argument
		,"PermutationSerial::permute(P_trans,x,y) : Error!" );
#endif
#ifdef ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY
	bool is_compatible;
	is_compatible = space_.is_compatible(x.space());
	THROW_EXCEPTION(
		!is_compatible, std::invalid_argument
		,"PermutationSerial::permute(P_trans,x,y) : Error, "
		"this->space().is_compatible(x.space()) returned false!" );
	is_compatible = space_.is_compatible(y->space());
	THROW_EXCEPTION(
		!is_compatible, std::invalid_argument
		,"PermutationSerial::permute(P_trans,x,y) : Error, "
		"this->space().is_compatible(y->space()) returned false!" );
#endif
	VectorDenseMutableEncap       y_d(*y);
	VectorDenseEncap              x_d(x);
	const IVector                 *p            = NULL;
	bool                          call_inv_perm = false;
	if( ( p = perm_.get() ) != NULL ) {
		if( P_trans == BLAS_Cpp::no_trans )
			call_inv_perm = false;
		else
			call_inv_perm = true;
	}
	else if( ( p = inv_perm_.get() ) != NULL ) {
		if( P_trans == BLAS_Cpp::no_trans )
			call_inv_perm = true;
		else
			call_inv_perm = false;
	}
	if( p ) {
		if( call_inv_perm )
			LinAlgPack::inv_perm_ele( x_d(), *p, &y_d() );
		else
			LinAlgPack::perm_ele( x_d(), *p, &y_d() );
	}
	else {
		// Just the identity permutation, nothing to do!
	}
}

void PermutationSerial::permute( 
	BLAS_Cpp::Transp          P_trans
	,VectorWithOpMutable      *y
	) const
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		y == NULL, std::invalid_argument
		,"PermutationSerial::permute(P_trans,y) : Error!" );
#endif
	VectorSpace::vec_mut_ptr_t
		t = y->clone();
	this->permute(P_trans,*t,y);
}

} // end namespace SparseLinAlgPack
