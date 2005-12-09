// //////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixConvertToSparseEncap.hpp
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

#ifndef MATRIX_CONVERT_TO_SPARSE_ENCAP_H
#define MATRIX_CONVERT_TO_SPARSE_ENCAP_H

#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace AbstractLinAlgPack {

///
/** Sparse conversion subclass based on views of a \c MatrixExtractSparseElements object.
 *
 * ToDo:Finish documentation!
 */
class MatrixConvertToSparseEncap
	: virtual public MatrixConvertToSparse
{
public:

	/** @name Public types */
	//@{

	///
	typedef Teuchos::RefCountPtr<const MatrixExtractSparseElements>  mese_ptr_t;
	///
	typedef Teuchos::RefCountPtr<const IVector>                      i_vector_ptr_t;

	//@}

	/** @name Constructors/initializers */
	//@{

	///
	/** Construct to uninitialized.
	 */
	MatrixConvertToSparseEncap();

	///
	/** Calls \c this->initialize().
	 */
	MatrixConvertToSparseEncap(
		const mese_ptr_t           &mese
		,const i_vector_ptr_t      &inv_row_perm
		,const Range1D             &row_rng
		,const i_vector_ptr_t      &inv_col_perm
		,const Range1D             &col_rng
		,const BLAS_Cpp::Transp    mese_trans
		,const value_type          alpha = 1.0
		);

	///
	/** Initialize a permuted view of a sparse matrix.
	 *
	 * <tt>A = alpha * op( (P'*B*Q)(row_rng,col_rng) )</tt>
	 *
	 * ToDo: Finish documentation!
	 */
	void initialize(
		const mese_ptr_t           &mese
		,const i_vector_ptr_t      &inv_row_perm
		,const Range1D             &row_rng
		,const i_vector_ptr_t      &inv_col_perm
		,const Range1D             &col_rng
		,const BLAS_Cpp::Transp    mese_trans
		,const value_type          alpha = 1.0
		);

	///
	/** Set uninitialized.
	 *
	 * ToDo: Finish documentation!
	 */
	void set_uninitialized();

	//@}

	/** @name Access */
	//@{

	///
	const mese_ptr_t& mese() const;
	///
	const i_vector_ptr_t& inv_row_perm() const;
	///
	const Range1D& row_rng() const;
	///
	const i_vector_ptr_t& inv_col_perm() const;
	///
	const Range1D& col_rng() const;
	///
	const BLAS_Cpp::Transp mese_trans() const;
	///
	const value_type alpha() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	size_type nz() const;

	//@}

	/** @name Overridden from MatrixConvertToSparse */
	//@{

	///
	index_type num_nonzeros(
		EExtractRegion        extract_region
		,EElementUniqueness   element_uniqueness
		) const;

	///
	void coor_extract_nonzeros(
		EExtractRegion                extract_region
		,EElementUniqueness           element_uniqueness
		,const index_type             len_Aval
		,value_type                   Aval[]
		,const index_type             len_Aij
		,index_type                   Arow[]
		,index_type                   Acol[]
		,const index_type             row_offset
		,const index_type             col_offset
		) const;

	//@}

private:

	typedef Teuchos::RefCountPtr<const VectorSpace> space_ptr_t;  

#ifdef DOXYGEN_COMPILE
	const MatrixExtractSparseElements    *mese;
	const DenseLinAlgPack::IVector            *inv_row_perm;
	Range1D                              row_rng;
	const DenseLinAlgPack::IVector            *inv_col_perm;
	Range1D                              col_rng;
#else
	mese_ptr_t          mese_;
	BLAS_Cpp::Transp    mese_trans_;
	value_type          alpha_;
	Range1D             row_rng_;
	Range1D             col_rng_;
	i_vector_ptr_t      inv_row_perm_;
	i_vector_ptr_t      inv_col_perm_;
	size_type           nz_full_;
	space_ptr_t         space_cols_;
	space_ptr_t         space_rows_;
#endif

};	// end class MatrixConvertToSparseEncap

// /////////////////////////////////////////////
// Inline members

// Access

inline
const MatrixConvertToSparseEncap::mese_ptr_t&
MatrixConvertToSparseEncap::mese() const
{
	return mese_;
}

inline
const MatrixConvertToSparseEncap::i_vector_ptr_t&
MatrixConvertToSparseEncap::inv_row_perm() const
{
	return inv_row_perm_;
}

inline
const Range1D& MatrixConvertToSparseEncap::row_rng() const
{
	return row_rng_;
}

inline
const MatrixConvertToSparseEncap::i_vector_ptr_t&
MatrixConvertToSparseEncap::inv_col_perm() const
{
	return inv_col_perm_;
}

inline
const Range1D& MatrixConvertToSparseEncap::col_rng() const
{
	return col_rng_;
}

inline
const BLAS_Cpp::Transp
MatrixConvertToSparseEncap::mese_trans() const
{
	return mese_trans_;
}

inline
const value_type MatrixConvertToSparseEncap::alpha() const
{
	return alpha_;
}
	
}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_CONVERT_TO_SPARSE_ENCAP_H
