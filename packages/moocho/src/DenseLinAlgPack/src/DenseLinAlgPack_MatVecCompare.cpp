// /////////////////////////////////////////////////////////////////////////////////
// MatVecCompare.cpp
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

#include <limits>
#include <algorithm>

#include "../include/MatVecCompare.h"
#include "../include/GenMatrixClass.h"
#include "../include/GenMatrixAsTriSym.h"

namespace {
using LinAlgPack::value_type;
using LinAlgPack::sqrt_eps;
inline
bool _comp(value_type val1, value_type val2)
{
	const value_type denom = std::_MAX( ::fabs(val1), 1.0 ); // compare relative errors.
	return ::fabs(val1 - val2) / denom < sqrt_eps;
} 

}	// end namespace

bool LinAlgPack::comp(const VectorSlice& vs1, const VectorSlice& vs2) {
	VectorSlice::const_iterator
		vs1_itr = vs1.begin(),
		vs2_itr = vs2.begin();
	for(; vs1_itr != vs1.end() && vs2_itr != vs2.end(); ++vs1_itr, ++vs2_itr)
		if( !_comp(*vs1_itr,*vs2_itr) ) return false;
	return true;
}

bool LinAlgPack::comp(const VectorSlice& vs, value_type alpha) {
	VectorSlice::const_iterator vs_itr = vs.begin();
	for(; vs_itr != vs.end(); ++vs_itr)
		if( !_comp(*vs_itr,alpha) ) return false;
	return true;
}

bool LinAlgPack::comp(const GenMatrixSlice& gms1, BLAS_Cpp::Transp trans1
	, const GenMatrixSlice& gms2, BLAS_Cpp::Transp trans2)
{
	for(size_type i = 1; i < std::_MIN(gms1.cols(),gms2.cols()); ++i)
		if( !comp( col(gms1,trans1,i) , col( gms2, trans2, i ) ) ) return false;
	return true;
}

bool LinAlgPack::comp(const GenMatrixSlice& gms, value_type alpha)
{
	for(size_type i = 1; i < gms.cols(); ++i)
		if( !comp( gms.col(i) , alpha ) ) return false;
	return true;
}

bool LinAlgPack::comp(const tri_ele_gms& tri_gms1, const tri_ele_gms& tri_gms2)
{
	using BLAS_Cpp::bool_to_trans;
	BLAS_Cpp::Transp
		trans1 = bool_to_trans(tri_gms1.uplo() != BLAS_Cpp::lower),
		trans2 = bool_to_trans(tri_gms2.uplo() != BLAS_Cpp::lower);

	const size_type n = tri_gms1.rows();	// same as cols()
	for(size_type i = 1; i < n; ++i) {
		if( !comp( col(tri_gms1.gms(),trans1,i)(i,n), col(tri_gms2.gms(),trans2,i)(i,n) ) )
			return false;
	}
	return true;
}

bool LinAlgPack::comp(const tri_ele_gms& tri_gms1, value_type alpha)
{
	using BLAS_Cpp::bool_to_trans;
	BLAS_Cpp::Transp
		trans1 = bool_to_trans(tri_gms1.uplo() != BLAS_Cpp::lower);

	const size_type n = tri_gms1.rows();	// same as cols()
	for(size_type i = 1; i < n; ++i) {
		if( !comp( col(tri_gms1.gms(),trans1,i)(i,n), alpha ) )
			return false;
	}
	return true;
}

bool LinAlgPack::comp_less(const VectorSlice& vs, value_type alpha)
{
	VectorSlice::const_iterator vs_itr = vs.begin();
	const value_type denom = std::_MAX( ::fabs(alpha), 1.0 );
	for(; vs_itr != vs.end(); ++vs_itr)
		if( *vs_itr > alpha ) return false;
	return true;
}
