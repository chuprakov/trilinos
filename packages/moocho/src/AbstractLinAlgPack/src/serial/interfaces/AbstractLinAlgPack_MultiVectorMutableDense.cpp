// /////////////////////////////////////////////////////////////////////////
// MultiVectorMutableDense.cpp
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

#include "SparseLinAlgPack/src/MultiVectorMutableDense.hpp"
#include "SparseLinAlgPack/src/VectorWithOpMutableDense.hpp"
#include "SparseLinAlgPack/src/MatrixSymWithOpGetGMSSymMutable.hpp"
#include "SparseLinAlgPack/src/SpVectorOp.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "DenseLinAlgPack/src/DMatrixOut.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "WorkspacePack.hpp"
#include "ThrowException.hpp"

namespace SparseLinAlgPack {

// Constructors / initializers

MultiVectorMutableDense::MultiVectorMutableDense(
	const size_type                    rows
	,const size_type                   cols
	)
{
	this->initialize(rows,cols);
}

MultiVectorMutableDense::MultiVectorMutableDense(
	DMatrixSlice                     gms
	,BLAS_Cpp::Transp                  gms_trans
	,const release_resource_ptr_t&     gms_release
	)
{
	this->initialize(gms,gms_trans,gms_release);
}

void MultiVectorMutableDense::initialize(
	const size_type                    rows
	,const size_type                   cols
	)
{
	namespace rcp = MemMngPack;
	namespace rmp = MemMngPack;
	typedef rcp::ref_count_ptr<DMatrix> vec_ptr_t;
	vec_ptr_t gms_ptr = rcp::rcp(new DMatrix(rows,cols));
	this->initialize(
		(*gms_ptr)()
		,BLAS_Cpp::no_trans
		,rcp::rcp(
			new rmp::ReleaseResource_ref_count_ptr<DMatrix>(
				gms_ptr
				)
			)
		);
}

void MultiVectorMutableDense::initialize(
	DMatrixSlice                     gms
	,BLAS_Cpp::Transp                  gms_trans
	,const release_resource_ptr_t&     gms_release
	)
{
	gms_.bind(gms);
	gms_trans_   = gms_trans;
	gms_release_ = gms_release;
}

// Overridden from MatrixWithOpGetGMS

const DMatrixSlice MultiVectorMutableDense::get_gms_view() const
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to create a copy and transpose it!
	}
	return get_gms(); // No memory to allocate!
}

void MultiVectorMutableDense::free_gms_view(const DMatrixSlice* gms_view) const
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to free the copy that we created in get_gms_view()
	}
	else {
		// Nothing to free!
	}
}

// Overridden from MatrixWithOpGetGMSMutable

DMatrixSlice MultiVectorMutableDense::get_gms_view()
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to create a copy and transpose it!
	}
	return set_gms(); // No memory to allocate!
}

void MultiVectorMutableDense::commit_gms_view(DMatrixSlice* gms_view)
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to free the copy that we created in get_gms_view()
	}
	else {
		// Nothing to free!
	}
}

// Overridden from MultiVector

MultiVectorMutableDense::access_by_t
MultiVectorMutableDense::access_by() const
{
	return ROW_ACCESS | COL_ACCESS | DIAG_ACCESS; // row, column and diagonal access is available!
}

// Overridden from MultiVectorMutable

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::col(index_type j)
{
	namespace rcp = MemMngPack;
	return rcp::rcp(
		new VectorWithOpMutableDense(
			DenseLinAlgPack::col( set_gms(), gms_trans(), j )
			,rcp::null ) );
}

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::row(index_type i)
{
	namespace rcp = MemMngPack;
	return rcp::rcp(
		new VectorWithOpMutableDense(
			DenseLinAlgPack::row( set_gms(), gms_trans(), i )
			,rcp::null ) );
}

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::diag(int k)
{
	namespace rcp = MemMngPack;
	return rcp::rcp(
		new VectorWithOpMutableDense(
			gms_.diag( gms_trans() == BLAS_Cpp::no_trans ? k : -k )
			,rcp::null ) );
}

MultiVectorMutableDense::multi_vec_mut_ptr_t
MultiVectorMutableDense::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
	namespace rcp = MemMngPack;
	return rcp::rcp(
		new MultiVectorMutableDense(
			gms_(
				gms_trans() == BLAS_Cpp::no_trans   ? row_rng : col_rng
				,gms_trans() == BLAS_Cpp::no_trans  ? col_rng : row_rng )
			,gms_trans()
			,rcp::null ) );
}

// Overridden from MatrixBase

size_type MultiVectorMutableDense::rows() const
{
	return BLAS_Cpp::rows( get_gms().rows(), get_gms().cols(), gms_trans() );
}

size_type MultiVectorMutableDense::cols() const
{
	return BLAS_Cpp::cols( get_gms().rows(), get_gms().cols(), gms_trans() );
}

// Overridden from MatrixOp

void MultiVectorMutableDense::zero_out()
{
	gms_ = 0.0;
}

void MultiVectorMutableDense::Mt_S( value_type alpha )
{
	DenseLinAlgPack::Mt_S(&gms_,alpha);
}

MatrixOp& MultiVectorMutableDense::operator=(const MatrixOp& mwo_rhs)
{
	DenseLinAlgPack::assign( &set_gms(), MatrixDenseEncap(mwo_rhs)(), gms_trans() );
	return *this;
}

std::ostream& MultiVectorMutableDense::output(std::ostream& out) const
{
	if(gms_trans() == BLAS_Cpp::no_trans)
		return out << gms_;
	return MatrixWithOpSerial::output(out);
}

bool MultiVectorMutableDense::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	,value_type beta, MatrixSymOp* sym_lhs
	) const
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	THROW_EXCEPTION(
		sym_lhs == NULL, std::invalid_argument
		,"MultiVectorMutableDense::syrk(...) : Error!" );
#endif
	MatrixSymWithOpGetGMSSymMutable
		*sym_get_lhs = dynamic_cast<MatrixSymWithOpGetGMSSymMutable*>(sym_lhs);
	if(!sym_get_lhs)
		return false;
	MatrixDenseSymMutableEncap  sym_gms_lhs(sym_get_lhs);
	DenseLinAlgPack::syrk( M_trans, alpha, get_gms(), beta, &sym_gms_lhs() );
	return true;
}

bool MultiVectorMutableDense::Mp_StMtM(
	MatrixOp* mwo_lhs, value_type alpha
	,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	,value_type beta ) const
{
	if(MultiVector::Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta))
		return true;
	return MatrixWithOpSerial::Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta);
}

bool MultiVectorMutableDense::Mp_StMtM(
	MatrixOp* mwo_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1
	,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	,value_type beta ) const
{
	if(MultiVector::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta))
		return true;
	return MatrixWithOpSerial::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
}

// Overridden from MatrixWithOpSerial

void MultiVectorMutableDense::Vp_StMtV(
	DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2, value_type beta) const
{
	DenseLinAlgPack::Vp_StMtV(
		vs_lhs,alpha,gms_,BLAS_Cpp::trans_trans(gms_trans(),trans_rhs1),vs_rhs2,beta);
}

void MultiVectorMutableDense::Vp_StMtV(
	DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	SparseLinAlgPack::Vp_StMtV(
		vs_lhs,alpha,gms_,BLAS_Cpp::trans_trans(gms_trans(),trans_rhs1),sv_rhs2,beta);
}

// protected

// Overridden from MultiVector

void MultiVectorMutableDense::apply_op(
	EApplyBy apply_by, const RTOpPack::RTOp& primary_op
	,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_objs[]
	,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
	,const index_type secondary_first_ele , const index_type secondary_sub_dim 
	) const
{
	MultiVector::apply_op(
		apply_by, primary_op, num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs
		,reduct_objs
		,primary_first_ele, primary_sub_dim, primary_global_offset, secondary_first_ele, secondary_sub_dim
		); // ToDo: Provide Specialized implementation if needed?
}

void MultiVectorMutableDense::apply_op(
	EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
	,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
	,const index_type secondary_first_ele , const index_type secondary_sub_dim 
	) const
{
	MultiVector::apply_op(
		apply_by, primary_op, secondary_op, num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs
		,reduct_obj
		,primary_first_ele, primary_sub_dim, primary_global_offset, secondary_first_ele, secondary_sub_dim
		); // ToDo: Provide Specialized implementation if needed?
}

} // end namespace SparseLinAlgPack
