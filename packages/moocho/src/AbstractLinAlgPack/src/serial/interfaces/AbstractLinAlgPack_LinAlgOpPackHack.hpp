// ///////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_LinAlgOpPackHack.hpp
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

#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

#include "DenseLinAlgPack/src/DenseLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp"

namespace LinAlgOpPack {

using DenseLinAlgPack::DVector;
using DenseLinAlgPack::DVectorSlice;
using DenseLinAlgPack::DMatrixSlice;
using AbstractLinAlgPack::SpVectorSlice;
using AbstractLinAlgPack::GenPermMatrixSlice;
using AbstractLinAlgPack::MatrixOp;
using AbstractLinAlgPack::MatrixNonsing;
using AbstractLinAlgPack::MatrixOpNonsing;

///
/** <tt>m_lhs += alpha * op(mwo_rhs1)</tt>.
 */
void Mp_StM(
	DMatrixSlice* vs_lhs, value_type alpha
	,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	);

///
/** <tt>vs_lhs = alpha * op(mwo_rhs1) * vs_rhs2 + beta * vs_lhs</tt>.
 */
void Vp_StMtV(
	DVectorSlice* vs_lhs, value_type alpha, const MatrixOp& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2
	,value_type beta = 1.0 );

///
/** <tt>vs_lhs = alpha * op(mwo_rhs1) * vs_rhs2 + beta * sv_lhs</tt>.
 */
void Vp_StMtV(
	DVectorSlice* vs_lhs, value_type alpha, const MatrixOp& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
	,value_type beta = 1.0 );

///
/** <tt>vs_lhs = inv(op(mwo_rhs1)) * vs_rhs2</tt>.
 */
void V_InvMtV(
	DVectorSlice* vs_lhs, const MatrixOpNonsing& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2 );

///
/** <tt>v_lhs = inv(op(mwo_rhs1)) * vs_rhs2</tt>.
 */
void V_InvMtV(
	DVector* v_lhs, const MatrixOpNonsing& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2 );

///
/** <tt>vs_lhs = inv(op(mwo_rhs1)) * sv_rhs2</tt>.
 */
void V_InvMtV(
	DVectorSlice* vs_lhs, const MatrixOpNonsing& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

///
/** <tt>v_lhs = inv(op(mwo_rhs1)) * sv_rhs2</tt>.
 */
void V_InvMtV(
	DVector* v_lhs, const MatrixOpNonsing& mwo_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

///
/** <tt>vs_lhs = alpha * op(gpms_rhs1) * op(mwo_rhs2) * vs_rhs3 + beta * vs_lhs</tt>.
 */
void Vp_StPtMtV(
	DVectorSlice* vs_lhs, value_type alpha
	,const GenPermMatrixSlice& gpms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	,const DVectorSlice& vs_rhs3, value_type beta = 1.0 );

///
/** <tt>vs_lhs = alpha * op(gpms_rhs1) * op(mwo_rhs2) * sv_rhs3 + beta * vs_lhs</tt>.
 */
void Vp_StPtMtV(
	DVectorSlice* vs_lhs, value_type alpha
	,const GenPermMatrixSlice& gpms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	,const SpVectorSlice& sv_rhs3, value_type beta = 1.0 );

} // end namespace LinAlgOpPack

#endif // LIN_ALG_OP_PACK_HACK_H
