// ///////////////////////////////////////////////////////////////////////
// LinAlgOpPackHack.hpp
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

#ifndef LIN_ALG_OP_PACK_HACK_H
#define LIN_ALG_OP_PACK_HACK_H

#include "SparseLinAlgPackTypes.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"

namespace LinAlgOpPack {

using LinAlgPack::Vector;
using LinAlgPack::VectorSlice;
using LinAlgPack::GenMatrixSlice;
using AbstractLinAlgPack::SpVectorSlice;
using AbstractLinAlgPack::GenPermMatrixSlice;

///
/** <tt>m_lhs += alpha * op(mwo_rhs1)</tt>.
 */
void Mp_StM(
	GenMatrixSlice* vs_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	);

///
/** <tt>vs_lhs = alpha * op(mwo_rhs1) * vs_rhs2 + beta * vs_lhs</tt>.
 */
void Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, const MatrixWithOp& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2
	,value_type beta = 1.0 );

///
/** <tt>vs_lhs = alpha * op(mwo_rhs1) * vs_rhs2 + beta * sv_lhs</tt>.
 */
void Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, const MatrixWithOp& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
	,value_type beta = 1.0 );

///
/** <tt>vs_lhs = inv(op(mwo_rhs1)) * vs_rhs2</tt>.
 */
void V_InvMtV(
	VectorSlice* vs_lhs, const MatrixWithOpNonsingular& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2 );

///
/** <tt>v_lhs = inv(op(mwo_rhs1)) * vs_rhs2</tt>.
 */
void V_InvMtV(
	Vector* v_lhs, const MatrixWithOpNonsingular& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2 );

///
/** <tt>vs_lhs = inv(op(mwo_rhs1)) * sv_rhs2</tt>.
 */
void V_InvMtV(
	VectorSlice* vs_lhs, const MatrixWithOpNonsingular& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

///
/** <tt>v_lhs = inv(op(mwo_rhs1)) * sv_rhs2</tt>.
 */
void V_InvMtV(
	Vector* v_lhs, const MatrixWithOpNonsingular& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

///
/** <tt>vs_lhs = alpha * op(gpms_rhs1) * op(mwo_rhs2) * vs_rhs3 + beta * vs_lhs</tt>.
 */
void Vp_StPtMtV(
	VectorSlice* vs_lhs, value_type alpha
	,const GenPermMatrixSlice& gpms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	,const VectorSlice& vs_rhs3, value_type beta = 1.0 );

///
/** <tt>vs_lhs = alpha * op(gpms_rhs1) * op(mwo_rhs2) * sv_rhs3 + beta * vs_lhs</tt>.
 */
void Vp_StPtMtV(
	VectorSlice* vs_lhs, value_type alpha
	,const GenPermMatrixSlice& gpms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	,const SpVectorSlice& sv_rhs3, value_type beta = 1.0 );

} // end namespace LinAlgOpPack

#endif // LIN_ALG_OP_PACK_HACK_H
