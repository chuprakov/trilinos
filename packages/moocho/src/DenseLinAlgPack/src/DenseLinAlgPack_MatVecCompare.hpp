// //////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_MatVecCompare.hpp
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

#ifndef MAT_VEC_COMPARE_H
#define MAT_VEC_COMPARE_H

#include <limits>
#if defined(_GNU_GXX)
#include <cmath>
#else
#include <math.h>
#endif

#include "DenseLinAlgPack_Types.hpp"
#include "TestingHelperPack_update_success.hpp"

namespace DenseLinAlgPack {

using TestingHelperPack::update_success;

/* * @name DVectorSlice and DMatrixSlice comparison functions.
  *
  * These functions compare the elements of two DVectorSlice or DMatrixSlice
  * objects.  If any of the corresponding elements does not obey
  * abs(ele1 - ele2) < sqrt(eps) then the functions return false, otherwise
  * they return true.  An exact test (bits) is not performed to allow for some round-off
  * error to occur and still equate.
  */
// @{

///
const value_type sqrt_eps
#if defined(_GNU_GXX)
	= std::sqrt(std::numeric_limits<value_type>::epsilon());
#elif defined(_CPQ_CXX)
	= ::sqrt(std::numeric_limits<value_type>::epsilon());
#else
	= ::sqrt(std::numeric_limits<value_type>::epsilon());
#endif

///
bool comp(const DVectorSlice& vs1, const DVectorSlice& vs2);

///
bool comp(const DVectorSlice& vs, value_type alpha);

///
bool comp(const DMatrixSlice& gms1, BLAS_Cpp::Transp trans1
	, const DMatrixSlice& gms2, BLAS_Cpp::Transp trans2);

/////
//bool comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2);

inline
///
bool comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2)
{
	return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
}

///
bool comp(const DMatrixSlice& gms1, value_type alpha);

///
bool comp(const DMatrixSliceTriEle& tri_gms1, const DMatrixSliceTriEle& tri_gms2);

///
bool comp(const DMatrixSliceTriEle& tri_gms1, value_type alpha);

///
bool comp_less(const DVectorSlice& vs, value_type alpha);

// @}

}	// end namespace DenseLinAlgPack

// ////////////////////////////////////
// Inline definitions

//inline
//bool DenseLinAlgPack::comp(const DMatrixSlice& gms1, const DMatrixSlice& gms2)
//{
//	return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
//}


#endif	// MAT_VEC_COMPARE_H
