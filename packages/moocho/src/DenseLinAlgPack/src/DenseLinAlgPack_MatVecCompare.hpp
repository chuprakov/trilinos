// //////////////////////////////////////////////////////////////////////
// MatVecCompare.h

#ifndef MAT_VEC_COMPARE_H
#define MAT_VEC_COMPARE_H

#include <limits>

#include "LinAlgPackTypes.h"
#include "Misc/include/update_success.h"

namespace LinAlgPack {

using TestingHelperPack::update_success;

/** @name VectorSlice and GenMatrixSlice comparison functions.
  *
  * These functions compare the elements of two VectorSlice or GenMatrixSlice
  * objects.  If any of the corresponding elements does not obey
  * abs(ele1 - ele2) < sqrt(eps) then the functions return false, otherwise
  * they return true.  An exact test (bits) is not performed to allow for some round-off
  * error to occure and still equate.
  */
//@{

///
const value_type sqrt_eps = ::sqrt(std::numeric_limits<value_type>::epsilon());

///
bool comp(const VectorSlice& vs1, const VectorSlice& vs2);

///
bool comp(const VectorSlice& vs, value_type alpha);

///
bool comp(const GenMatrixSlice& gms1, BLAS_Cpp::Transp trans1
	, const GenMatrixSlice& gms2, BLAS_Cpp::Transp trans2);

/////
//bool comp(const GenMatrixSlice& gms1, const GenMatrixSlice& gms2);

inline
///
bool comp(const GenMatrixSlice& gms1, const GenMatrixSlice& gms2)
{
	return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
}

///
bool comp(const GenMatrixSlice& gms1, value_type alpha);

///
bool comp(const tri_ele_gms& tri_gms1, const tri_ele_gms& tri_gms2);

///
bool comp(const tri_ele_gms& tri_gms1, value_type alpha);

///
bool comp_less(const VectorSlice& vs, value_type alpha);

//@}

}	// end namespace LinAlgPack

// ////////////////////////////////////
// Inline definitions

//inline
//bool LinAlgPack::comp(const GenMatrixSlice& gms1, const GenMatrixSlice& gms2)
//{
//	return comp(gms1, BLAS_Cpp::no_trans, gms2, BLAS_Cpp::no_trans);
//}


#endif	// MAT_VEC_COMPARE_H