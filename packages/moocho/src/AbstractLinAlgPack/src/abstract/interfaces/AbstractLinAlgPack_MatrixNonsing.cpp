// //////////////////////////////////////////////////////////////////////////////////
//  MatrixNonsingular.cpp
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

// ToDo: 3/6/00: Provide default implementations for these
// operations.

#include <assert.h>

#include "AbstractLinAlgPack/include/MatrixNonsingular.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/SpVectorView.h"
#include "AbstractLinAlgPack/include/EtaVector.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

namespace AbstractLinAlgPack {

// Clone

MatrixNonsingular::mat_mns_mut_ptr_t
MatrixNonsingular::clone_mns()
{
	return MemMngPack::null;
}

MatrixNonsingular::mat_mns_ptr_t
MatrixNonsingular::clone_mns() const
{
	return const_cast<MatrixNonsingular*>(this)->clone_mns(); // Implicit conversion to const
}

// Level-2 BLAS

void MatrixNonsingular::V_InvMtV(
	VectorWithOpMutable* y, BLAS_Cpp::Transp M_trans, const SpVectorSlice& sx
	) const
{
	if( sx.nz() ) {
		VectorSpace::vec_mut_ptr_t
			x = (M_trans == BLAS_Cpp::no_trans
					  ? this->space_cols()
					  : this->space_rows()
				).create_member();
		x->set_sub_vector(sub_vec_view(sx));
		this->V_InvMtV(y,M_trans,*x);
	}
	else {
		*y = 0.0;
	}
}

value_type MatrixNonsingular::transVtInvMtV(
	const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3
	) const
{
	VectorSpace::vec_mut_ptr_t
		v = (trans_rhs2 == BLAS_Cpp::no_trans
				  ? this->space_rows()
				  : this->space_cols()
			).create_member();
	this->V_InvMtV( v.get(), trans_rhs2, v_rhs3 );
	return dot(v_rhs1,*v);
}

value_type MatrixNonsingular::transVtInvMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3
	) const
{
	VectorSpace::vec_mut_ptr_t
		v = (trans_rhs2 == BLAS_Cpp::no_trans
				  ? this->space_rows()
				  : this->space_cols()
			).create_member();
	this->V_InvMtV( v.get(), trans_rhs2, sv_rhs3 );
	return dot(sv_rhs1,*v);
}

// Level-3 BLAS

void MatrixNonsingular::M_StInvMtM(
	MatrixWithOp* C_lhs, value_type alpha
	,BLAS_Cpp::Transp M_trans
	,const MatrixWithOp& B, BLAS_Cpp::Transp B_trans
	) const
{
	//
	// C = a * inv(op(M)) * op(B)
	//
	using DynamicCastHelperPack::dyn_cast;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
#ifdef _DEBUG
	THROW_EXCEPTION(
		C_lhs == NULL, std::invalid_argument
		,"MatrixNonsingular::M_StInvMtM(...) : Error!" );
	
#endif
	const size_type
		C_rows = C_lhs->rows(),
		C_cols = C_lhs->cols();
	const size_type
		op_B_cols = BLAS_Cpp::cols( B.rows(), B.cols(), B_trans );
#ifdef _DEBUG
	// We can't check vector spaces since *this may not support MatrixWithOp
	// However, we could dynamic cast to see if MatrixWithOp is supported and then
	// be able to use Mp_MtM_assert_compatibility() but this is okay for now.
	const size_type
		M_rows    = this->rows(),
		M_cols    = this->cols(),
		op_B_rows = BLAS_Cpp::rows( B.rows(), B.cols(), B_trans );
	THROW_EXCEPTION(
		C_rows != M_rows || M_rows != M_cols || M_cols != op_B_rows || C_cols != op_B_cols
		, std::invalid_argument
		,"MatrixNonsingular::M_StInvMtM(...) : Error!" );
#endif
	//
	// Compute C = a * inv(op(M)) * op(B) one column at a time:
	//
	// C(:,j) = inv(op(M)) * a * op(B) * e(j)    , for j = 1...C.cols()
	//                       \______________/    
	//                              t_j
	//
	MultiVectorMutable  &C = dyn_cast<MultiVectorMutable>(*C_lhs);
	VectorSpace::vec_mut_ptr_t
		t_j = ( B_trans == no_trans ? B.space_cols() : B.space_rows() ).create_member();
	for( size_type j = 1; j <= C_cols; ++j ) {
		// t_j = alpha * op(B) * e_j
		EtaVector e_j( j, op_B_cols );
		LinAlgOpPack::V_StMtV( t_j.get(), alpha, B, B_trans, e_j() );
		// C(:,j) = inv(op(M)) * t_j
		AbstractLinAlgPack::V_InvMtV( C.col(j).get(), *this, M_trans, *t_j );
	}
}

void MatrixNonsingular::M_StMtInvM(
	MatrixWithOp* g_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,BLAS_Cpp::Transp trans_rhs2
	) const
{
	assert(0); // ToDo: Implement!
}

}	// end namespace AbstractLinAlgPack
