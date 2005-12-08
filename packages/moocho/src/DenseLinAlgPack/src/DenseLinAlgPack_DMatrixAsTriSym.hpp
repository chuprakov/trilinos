// //////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_DMatrixAsTriSym.hpp
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

#ifndef GEN_MATRIX_AS_TRI_SYM_H
#define GEN_MATRIX_AS_TRI_SYM_H

#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

/* * @name	Packaging arguments for GenMatrixSlices treated as triangular
  *			and symmetric matrices in BLAS-like linear algebra operations.
  *
  */

// @{

// /////////////////////////////////////////////////////////////////////////////////////
///
/* * Aggregate information for a triangular matrix (element-wise) stored in a DMatrix.
  *
  * This is the type to be used as lhs and rhs arguments in element-wise
  * linear algebra operations like assignment and binary arithmetic.
  */
class DMatrixSliceTriEle {
public:
	///
	DMatrixSliceTriEle(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
		: gms_(const_cast<DMatrixSlice&>(gms)), uplo_(uplo)
	{
		#ifdef LINALGPACK_CHECK_RHS_SIZES
			assert_gms_square(gms);
		#endif
	}
	///
	size_type rows() const {
		return gms_.rows();
	}
	///
	size_type cols() const {
		return gms_.cols();
	}
	///
	DMatrixSlice& gms() {
		return gms_;
	}
	///
	const DMatrixSlice& gms() const {
		return gms_;
	}
	///
	BLAS_Cpp::Uplo	uplo() const {
		return uplo_;
	}
	/// Allow address to be taken of rvalue of this object
	DMatrixSliceTriEle* operator&() {
	  return this;
	}
	///
	const DMatrixSliceTriEle* operator&() const {
	  return this;
	}

private:	
	DMatrixSlice	gms_;
	BLAS_Cpp::Uplo	uplo_;
	// Not defined and not to be called
	DMatrixSliceTriEle();
	DMatrixSliceTriEle& operator=(const DMatrixSliceTriEle&);
};	// end class DMatrixSliceTriEle

inline
/// Return a triangular element-wise matrix
DMatrixSliceTriEle nonconst_tri_ele(DMatrixSlice gms, BLAS_Cpp::Uplo uplo)
{
	return DMatrixSliceTriEle(gms, uplo);
}

inline
///
const DMatrixSliceTriEle tri_ele(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
	return DMatrixSliceTriEle(gms, uplo);
}

// ////////////////////////////////////////////////////////////////////////////////
///
/* * Aggregate information for a triangular matrix (structure dependent) stored in a DMatrix.
  *
  * This is the type to be used as a rhs argument in linear algebra operations
  * that are structure specific like the BLAS operations.
  */
class DMatrixSliceTri {
public:
	///
	DMatrixSliceTri(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
		: gms_(const_cast<DMatrixSlice&>(gms)), uplo_(uplo), diag_(diag)
	{
		#ifdef LINALGPACK_CHECK_RHS_SIZES
			assert_gms_square(gms);
		#endif
	}
	///
	size_type rows() const {
		return gms_.rows();
	}
	///
	size_type cols() const {
		return gms_.cols();
	}
	///
	DMatrixSlice& gms() {
		return gms_;
	}
	///
	const DMatrixSlice& gms() const {
		return gms_;
	}
	///
	BLAS_Cpp::Uplo	uplo() const {
		return uplo_;
	}
	///
	BLAS_Cpp::Diag	diag() const {
		return diag_;
	}
	/// Allow address to be taken of rvalue of this object
	DMatrixSliceTri* operator&() {
	  return this;
	}
	///
	const DMatrixSliceTri* operator&() const {
	  return this;
	}

private:	
	DMatrixSlice	gms_;
	BLAS_Cpp::Uplo	uplo_;
	BLAS_Cpp::Diag	diag_;
	// not defined and not to be called
	DMatrixSliceTri();
	DMatrixSliceTri& operator=(const DMatrixSliceTri&);
};	// end class DMatrixSliceTri

inline
/// Return a triangular matrix
DMatrixSliceTri nonconst_tri(DMatrixSlice gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
{
	return DMatrixSliceTri(gms, uplo, diag);
}

inline
///
const DMatrixSliceTri tri(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
{
	return DMatrixSliceTri(gms, uplo, diag);
}

// /////////////////////////////////////////////////////////////////////////////////////////
///
/* * Aggregate information for a symmetric matrix stored in a DMatrix.
  *
  * This is the type to be used as both lhs and rhs arguments in linear algebra operations.
  */
class DMatrixSliceSym {
public:
	///
	DMatrixSliceSym(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
		: gms_(const_cast<DMatrixSlice&>(gms)), uplo_(uplo)
	{
		#ifdef LINALGPACK_CHECK_RHS_SIZES
			assert_gms_square(gms);
		#endif
	}
	///
	size_type rows() const {
		return gms_.rows();
	}
	///
	size_type cols() const {
		return gms_.cols();
	}
	///
	DMatrixSlice& gms() {
		return gms_;
	}
	///
	const DMatrixSlice& gms() const {
		return gms_;
	}
	///
	BLAS_Cpp::Uplo	uplo() const {
		return uplo_;
	}
	/// Allow address to be taken of rvalue of this object
	DMatrixSliceSym* operator&() {
	  return this;
	}
	const DMatrixSliceSym* operator&() const {
	  return this;
	}

private:	
	DMatrixSlice	gms_;
	BLAS_Cpp::Uplo	uplo_;
	// not defined and not to be called
	DMatrixSliceSym();
	DMatrixSliceSym& operator=(const DMatrixSliceTri&);
};	// end class DMatrixSliceSym

inline
/// Return a symmetric matrix
DMatrixSliceSym nonconst_sym(DMatrixSlice gms, BLAS_Cpp::Uplo uplo)
{
	return DMatrixSliceSym(gms, uplo);
}

inline
///
const DMatrixSliceSym sym(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
	return DMatrixSliceSym(gms, uplo);
}

// @}

} // end namespace DenseLinAlgPack

#endif	// GEN_MATRIX_AS_TRI_SYM_H
