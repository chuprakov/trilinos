// ///////////////////////////////////////////////////////////////////////
// LinAlgOpPackHack.cpp
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
//

#include "SparseLinAlgPack/include/LinAlgOpPackHack.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"

void LinAlgOpPack::Vp_StMtV(
	VectorSlice* y, value_type a, const MatrixWithOp& M
	,BLAS_Cpp::Transp M_trans, const VectorSlice& x, value_type b )
{
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::VectorDenseMutableEncap;
	VectorSpace::vec_mut_ptr_t
		ay = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member(),
		ax = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member();
	(VectorDenseMutableEncap(*ay))() = *y;
	(VectorDenseMutableEncap(*ax))() = x;
	Vp_StMtV( ay.get(), a, M, M_trans, *ax, b );
	*y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::Vp_StMtV(
	VectorSlice* y, value_type a, const MatrixWithOp& M
	,BLAS_Cpp::Transp M_trans, const SpVectorSlice& x, value_type b )
{
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::VectorDenseMutableEncap;
	VectorSpace::vec_mut_ptr_t
		ay = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member();
	(VectorDenseMutableEncap(*ay))() = *y;
	Vp_StMtV( ay.get(), a, M, M_trans, x, b );
	*y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
	VectorSlice* y, const MatrixWithOpNonsingular& M
	,BLAS_Cpp::Transp M_trans, const VectorSlice& x )
{
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::VectorDenseMutableEncap;
	VectorSpace::vec_mut_ptr_t
		ay = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member(),
		ax = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member();
	(VectorDenseMutableEncap(*ay))() = *y;
	(VectorDenseMutableEncap(*ax))() = x;
	V_InvMtV( ay.get(), M, M_trans, *ax );
	*y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
	Vector* y, const MatrixWithOpNonsingular& M
	,BLAS_Cpp::Transp M_trans, const VectorSlice& x )
{
	y->resize( BLAS_Cpp::rows( M.rows(), M.cols(), M_trans ) );
	V_InvMtV( &(*y)(), M, M_trans, x );
}

void LinAlgOpPack::V_InvMtV(
	VectorSlice* y, const MatrixWithOpNonsingular& M
	,BLAS_Cpp::Transp M_trans, const SpVectorSlice& x )
{
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::VectorDenseMutableEncap;
	VectorSpace::vec_mut_ptr_t
		ay = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member();
	(VectorDenseMutableEncap(*ay))() = *y;
	V_InvMtV( ay.get(), M, M_trans, x );
	*y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
	Vector* y, const MatrixWithOpNonsingular& M
	,BLAS_Cpp::Transp M_trans, const SpVectorSlice& x )
{
	y->resize( BLAS_Cpp::rows( M.rows(), M.cols(), M_trans ) );
	V_InvMtV( &(*y)(), M, M_trans, x );
}

void LinAlgOpPack::Vp_StPtMtV(
	VectorSlice* y, value_type a
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,const MatrixWithOp& M, BLAS_Cpp::Transp M_trans
	,const VectorSlice& x, value_type b )
{
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::VectorDenseMutableEncap;
	using AbstractLinAlgPack::Vp_StPtMtV;
	VectorSpace::vec_mut_ptr_t
		ay = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member(),
		ax = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member();
	(VectorDenseMutableEncap(*ay))() = *y;
	(VectorDenseMutableEncap(*ax))() = x;
	Vp_StPtMtV( ay.get(), a, P, P_trans, M, M_trans, *ax, b );
	*y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::Vp_StPtMtV(
	VectorSlice* y, value_type a
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,const MatrixWithOp& M, BLAS_Cpp::Transp M_trans
	,const SpVectorSlice& x, value_type b )
{
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::VectorDenseMutableEncap;
	using AbstractLinAlgPack::Vp_StPtMtV;
	VectorSpace::vec_mut_ptr_t
		ay = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member();
	(VectorDenseMutableEncap(*ay))() = *y;
	Vp_StPtMtV( ay.get(), a, P, P_trans, M, M_trans, x, b );
	*y = VectorDenseMutableEncap(*ay)();
}
