// /////////////////////////////////////////////////////////////////////////
// MultiVectorMutableDense.h
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

#ifndef MULTI_VECTOR_MUTABLE_DENSE_H
#define MULTI_VECTOR_MUTABLE_DENSE_H

#include "MatrixWithOpSerial.h"
#include "MatrixWithOpGetGMSMutable.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "ReleaseResource.h"

namespace SparseLinAlgPack {

///
/** <tt>MultiVectorMutable</tt> "Adapter" subclass for <tt>LinAlgPack::GenMatrixSlice</tt>
 * or <tt>LinAlgPack::GenMatrix</tt> object.
 *
 * This class can be used either as a view of a <tt>LinAlgPack::GenMatrixSlice</tt> object
 * or as a storage type for a <tt>LinAlgPack::GenMatrix</tt> object.
 *
 * To create a storage type with the dimensions of <tt>rows x cols</tt> just call the
 * constructor <tt>MatrixWithOpMutableDense(rows,cols)</tt> or after construction you
 * can call <tt>this->initialize(rows,cols)</tt>.
 *
 * To simply create a view of a matrix or its transpose, say \c M, without ownership
 * just call <tt>MatrixWithOpMutableDense(M(),M_trans,NULL)</tt> or after construction call
 * <tt>this->initialize(M(),M_trans,NULL)</tt>.
 *
 * Alternately, \c this can be given a matrix with the responsibility to delete any associated
 * memory by calling <tt>this->initialize()</tt> with a <tt>ReleaseResource</tt> object to
 * perform the deallocation.
 *
 * If \c this has been initialized by <tt>this->initialize(rows,cols)</tt> and
 * if the client really needs to get at the <tt>LinAlgPack::GenMatrix</tt> object
 * itself, then it can be obtained as:
 \code
 void f( MatrixWithOpMutableDense* M )
     namespace rmp = MemMngPack;
     GenMatrix &_M = *dynamic_cast<rmp::ReleaseResource_ref_count_ptr<GenMatrix&> >(*M.gms_release()).ptr;
 \endcode
 * This is not pretty but it is not supposed to be.  Of course the above function will throw
 * an exception if the <tt>dynamic_cast<></tt> fails.
 */
class MultiVectorMutableDense
	: public AbstractLinAlgPack::MultiVectorMutable   // doxygen needs the full path
	, public MatrixWithOpSerial
	, public MatrixWithOpGetGMS
{
public:

	///
	typedef MemMngPack::ref_count_ptr<
		MemMngPack::ReleaseResource>  release_resource_ptr_t;

	/** @name Constructors / initializers */
	//@{

	///
	/** Calls <tt>this->initialize(rows,cols)</tt>.
	 */
	MultiVectorMutableDense(
		const size_type                    rows = 0
		,const size_type                   cols = 0
		);
	///
	/** Calls <tt>this->initialize(gms,gms_trans,gms_release)</tt>.
	 */
	MultiVectorMutableDense(
		GenMatrixSlice                     gms
		,BLAS_Cpp::Transp                  gms_trans
		,const release_resource_ptr_t&     gms_release
		);
	///
	/** Call <tt>this->initialize(v,v_release)</tt> with an allocated <tt>LinAlgPack::Vector</tt>
	 * object.
	 */
	void initialize(
		const size_type                    rows
		,const size_type                   cols
		);
	///
	/** Initialize with a dense matrix slice.
	 *
	 * Note that solve of the method overrides have to allocate a temporary memory
	 * if <tt>gms_trans != no_trans</tt> (see \c get_gms_view() and \c output()).
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->rows() = ( gms_trans = no_trans ? gms.rows() : gms.cols())</tt>
	 * <li> <tt>this->cols() = ( gms_trans = no_trans ? gms.cols() : gms.rows())</tt>
	 * </ul>
	 */
	void initialize(
		GenMatrixSlice                     gms
		,BLAS_Cpp::Transp                  gms_trans
		,const release_resource_ptr_t&     gms_release
		);

	//@}

	/** @name Access */
	//@{
	
	///
	/** Return the non-const dense matrix.
	 *
	 * Note that calling this method may result in the matrix implementation
	 * being modified.  Therefore, no other methods on \c this object should be
	 * called until the <tt>GenMatrixSlice</tt> returned from this method is
	 * discarded.
	 */
	GenMatrixSlice set_gms();
	///
	/** Return a const dense matrix.
	 */
	const GenMatrixSlice get_gms() const;
	///
	/** Return if underlying matrix is being viewed as the transpose or non-transposed.
	 */
	BLAS_Cpp::Transp gms_trans() const;
	///
	/** Return a <tt>ref_count_ptr<></tt> pointer to the object that will
	 * release the associated resource.
	 */
	const release_resource_ptr_t& gms_release() const;

	//@}

	/** @name Overridden from MatrixWithOpGetGMS */
	//@{

	///
	const GenMatrixSlice get_gms_view() const;
	///
	void free_gms_view(const GenMatrixSlice* gms_view) const;

	//@}

	/** @name Overridden from MatrixWithOpGetGMSMutable */
	//@{

	///
	GenMatrixSlice get_gms_view();
	///
	void commit_gms_view(GenMatrixSlice* gms_view);

	//@}
	
	/** @name Overridden from MatrixBase */
	//@{

	///
	size_type rows() const;
	///
	size_type cols() const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	void zero_out();
	///
	void Mt_S( value_type alpha );
	///
	MatrixWithOp& operator=(const MatrixWithOp& mwo_rhs);
	///
	std::ostream& output(std::ostream& out) const;

	//@}

	/** @name Overridden from MultiVector */
	//@{

	///
	access_by_t access_by() const;

	//@}

	/** @name Overridden from MultiVectorMutable */
	//@{

	///
	vec_mut_ptr_t col(index_type j);
	///
	vec_mut_ptr_t row(index_type i);
	///
	vec_mut_ptr_t diag(int k);
	///
	multi_vec_mut_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng);

	//@}

	/** @name Overridden from MatrixWithOpSerial */
	//@{

	///
	void Vp_StMtV(
		VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	// ToDo: Add more overrides as they are needed!

	//@}

protected:

	/** @name Overridden from MultiVector */
	//@{

	///
	void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_objs[]
		,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
		,const index_type secondary_first_ele , const index_type secondary_sub_dim 
		) const;
	///
	void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
		,const index_type secondary_first_ele , const index_type secondary_sub_dim 
		) const;

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	GenMatrixSlice            gms_;
	BLAS_Cpp::Transp          gms_trans_;
	release_resource_ptr_t    gms_release_;

}; // end class MultiVectorMutableDense			

// //////////////////////////////////////
// Inline members

inline
GenMatrixSlice
MultiVectorMutableDense::set_gms()
{
	return gms_;
}

inline
const GenMatrixSlice
MultiVectorMutableDense::get_gms() const
{
	return gms_;
}

inline
BLAS_Cpp::Transp
MultiVectorMutableDense::gms_trans() const
{
	return gms_trans_;
}

inline
const MultiVectorMutableDense::release_resource_ptr_t&
MultiVectorMutableDense::gms_release() const
{
	return gms_release_;
}

} // end namespace SparseLinAlgPack

#endif // MULTI_VECTOR_MUTABLE_DENSE_H
