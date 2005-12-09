// ///////////////////////////////////////////////////////
// ConstrainedOptPack_MatrixGenBanded.hpp
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

#ifndef MATRIX_GEN_BANDED_H
#define MATRIX_GEN_BANDED_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Miref_count_ptr.h"
#include "MiReleaseResource.h"

namespace ConstrainedOptPack {
///
/** Matrix subclass for general (possibly singular) banded matrices.
 * 
 * The banded matrix is stored by column in a simple flat rectangular matrix.
 * For example, for #m = 10, n = 8, kl = 3, ku = 2# the matrix #M# is stored in the
 * following format #MB# (same as for the BLAS routine xGBMV(...)):
 \begin{verbatim}

         M                                   MB
 [ x x x               ]
 [ x x x x             ]         [ o o x x x x x x x x ] \ ku = 2
 [ x x x x x           ]         [ o x x x x x x x x o ] /
 [ x x x x x x         ]    =>   [ x x x x x x x x o o ]
 [   x x x x x x       ]         [ x x x x x x x o o o ] \
 [     x x x x x x     ]         [ x x x x x x o o o o ] | kl = 3
 [       x x x x x x   ]         [ x x x x x o o o o o ] /
 [         x x x x x x ]           1 2 3 4 5 6 7 8 9 0
   1 2 3 4 5 6 7 8 9 0

 \end{verbatim}
 */
class MatrixGenBanded	: public MatrixOp
{
public:
	
	///
	typedef Teuchos::RefCountPtr<
		MemMngPack::ReleaseResource>  release_resource_ptr_t;

	// //////////////
    // Constructors

	///
	/** Construct and Initialize.
	 *
	 * This constructor just calls #this->initialize(...)#.
	 */
	MatrixGenBanded(
		size_type                         m                       = 0
		,size_type                        n                       = 0
		,size_type                        kl                      = 0
		,size_type                        ku                      = 0
		,DMatrixSlice                   *MB                     = NULL
		,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
		);

	// ///////////////////////////
    // Access representation

	///
	/** Initialize
	 *
	 * If called with all of the default arguments then #this# will become uninitialized.
	 *
	 * ToDo: Finish pre and post conditions!
	 *
	 * @param  m        [in] Determines the size of the banded matrix (m x n).  If
	 *                  If #m == 0# then all of the following arguments should be left at
	 * @param  n        [in] Determines the size of the banded matrix (m x n).
	 * @param  kl       [in] Determines the lower band width of the matrix as defined by xGBMV(...).
	 * @param  ku       [in] Determines the band width of the matrix as defined by xGBMV(...).
	 * @param  MB       [in/state] If #MB != NULL# then this matrix (size (kl+ku+1) x n) is used to store
	 *                  the original banded matrix #M# in the format of xGBMV(...).  This matrix must
	 *                  be initialized on input.
	 * @param  MB_release_resource_ptr
	 *                  [in] Only significant if #MB != NULL#.  Points to a resource to
	 *                  be released when #MB# is no longer needed.
	 */
	void initialize(
		size_type                         m                       = 0
		,size_type                        n                       = 0
		,size_type                        kl                      = 0
		,size_type                        ku                      = 0
		,DMatrixSlice                   *MB                     = NULL
		,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
		);

	///
	size_type kl() const;
	///
	size_type ku() const;
	///
	/** Get view of MB.
	 */
	DMatrixSlice& MB();
	///
	const DMatrixSlice& MB() const;

	// /////////////////////////////
	// Overridden from MatrixOp

	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	size_type nz() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const DVectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const DVectorSlice& vs_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;

private:
	
	// /////////////////////////////
	// Private data members

	size_type                       m_;
	size_type                       n_;
	size_type                       kl_;
	size_type                       ku_;
	DMatrixSlice                  MB_;
	release_resource_ptr_t          MB_release_resource_ptr_;

	// /////////////////////////////
	// Private member functions

	void assert_initialized() const;

}; // end class MatrixGenBanded

// ///////////////////////////////////////////////////////
// Inline members for MatrixGenBanded

inline
size_type MatrixGenBanded::kl() const
{
	return kl_;
}

inline
size_type MatrixGenBanded::ku() const
{
	return ku_;
}

inline
DMatrixSlice& MatrixGenBanded::MB()
{
	return MB_;
}

inline
const DMatrixSlice& MatrixGenBanded::MB() const
{
	return MB_;
}

} // end namespace ConstrainedOptPack

#endif // MATRIX_GEN_BANDED_H
