// //////////////////////////////////////////////////////////////////////////////////
// LinAlgOpPackDef.h
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

#ifndef ABSTRACT_LIN_ALG_OP_PACK_DEF_H
#define ABSTRACT_LIN_ALG_OP_PACK_DEF_H

#include "LinAlgOpPackDecl.h"	// also includes some inline function definitions
#include "MatrixWithOp.h"
#include "VectorWithOpMutable.h"
#include "VectorStdOps.h"
#include "LinAlgPackAssertOp.h"

namespace LinAlgOpPack {

using BLAS_Cpp::rows;
using BLAS_Cpp::cols;

// Inject assert functions
using AbstractLinAlgPack::Vp_V_assert_sizes;
using AbstractLinAlgPack::VopV_assert_sizes;
using AbstractLinAlgPack::Mp_M_assert_sizes;
using AbstractLinAlgPack::MopM_assert_sizes;
using AbstractLinAlgPack::Vp_MtV_assert_sizes;
using AbstractLinAlgPack::MtV_assert_sizes;
using AbstractLinAlgPack::MtM_assert_sizes;

// Inject names of base linear algebra functions for LinAlgPack.
// Note that this is neccesary in MS VC++ 5.0 because
// it does not perform name lookups properly but it
// is not adverse to the standard so it is a portable
// fix.
//using AbstractLinAlgPack::Vt_S;
using AbstractLinAlgPack::Vp_StV;
using AbstractLinAlgPack::Vp_StMtV;
//using AbstractLinAlgPack::Mt_S;
using AbstractLinAlgPack::Mp_StM;
using AbstractLinAlgPack::Mp_StMtM;

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Vectors

// v_lhs = V_rhs.
template <class V>
void assign(VectorWithOpMutable* v_lhs, const V& V_rhs) {
	Vp_V_assert_sizes( v_lhs->dim(), V_rhs.dim() );
	(*v_lhs) = 0.0;
	Vp_V(v_lhs,V_rhs);
}

// v_lhs = alpha * V_rhs.
template <class V>
void V_StV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs) {
	Vp_V_assert_sizes( v_lhs->dim(), V_rhs.dim() );
	(*v_lhs) = 0.0;
	Vp_StV(v_lhs,alpha,V_rhs);
}

// v_lhs = V1_rhs1 + V2_rhs2.
template <class V1, class V2>
void V_VpV(VectorWithOpMutable* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
	Vp_V_assert_sizes( v_lhs->dim(), V1_rhs1.dim() );
	(*v_lhs) = 0.0;
	Vp_V(v_lhs,V1_rhs1);
	Vp_V(v_lhs,V2_rhs2);
}

// v_lhs = V_rhs1 - V_rhs2.
template <class V1, class V2>
void V_VmV(VectorWithOpMutable* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
	Vp_V_assert_sizes( v_lhs->dim(), V1_rhs1.dim() );
	(*v_lhs) = 0.0;
	Vp_V(v_lhs,V1_rhs1);
	Vp_StV(v_lhs,-1.0,V2_rhs2);
}

// v_lhs = alpha * V_rhs1 + v_rhs2.
template <class V>
void V_StVpV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs1
	, const VectorWithOp& v_rhs2)
{
	VopV_assert_sizes(V_rhs1.dim(),v_rhs2.dim());
	(*v_lhs) = v_rhs2;
	Vp_StV(v_lhs,alpha,V_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Matrices

// M_lhs = op(M_rhs).
template <class M>
void assign(MatrixWithOp* M_lhs, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	Mp_M_assert_sizes(M_lhs->rows(), M_lhs->cols(), BLAS_Cpp::no_trans
		, M_rhs.rows(), M_rhs.cols(), trans_rhs	);
	M_lhs->zero_out();
	Mp_StM(M_lhs,1.0,M_rhs,trans_rhs);
}

// M_lhs = alpha * op(M_rhs).
template <class M>
void M_StM(MatrixWithOp* M_lhs, value_type alpha, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	Mp_M_assert_sizes(M_lhs->rows(), M_lhs->cols(), BLAS_Cpp::no_trans
		, M_rhs.rows(), M_rhs.cols(), trans_rhs	);
	M_lhs->zero_out();
	Mp_StM(M_lhs,alpha,M_rhs,trans_rhs);
}

// M_lhs = op(M1_rhs1) + op(M2_rhs2).
template <class M1, class M2>
void M_MpM(MatrixWithOp* M_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
						,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
	assert_M_lhs(*M_lhs, rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
						   , cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
	M_lhs->zero_out();
	Mp_M(M_lhs,M1_rhs1,trans_rhs1);
	Mp_M(M_lhs,M2_rhs2,trans_rhs2);
}

// M_lhs = op(M_rhs1) - op(M_rhs2).
template <class M1, class M2>
void M_MmM(MatrixWithOp* M_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
						,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
	assert_M_lhs(*M_lhs, rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
						   , cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
	M_lhs->zero_out();
	Mp_M(M_lhs,M1_rhs1,trans_rhs1);
	Mp_StM(M_lhs,-1.0,M2_rhs2,trans_rhs2);
}

// M_lhs = alpha * op(M_rhs1) + op(m_rhs2).
template <class M>
void M_StMpM(MatrixWithOp* M_lhs, value_type alpha, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M_rhs1.rows(),M_rhs1.cols(),trans_rhs1
						,m_rhs2.rows(),m_rhs2.cols(),trans_rhs2);
	assign(M_lhs,m_rhs2,trans_rhs2);
	Mp_StM(M_lhs,alpha,M_rhs1,trans_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////// /////
// Level 2 BLAS

// v_lhs = alpha * op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, const M& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2)
{
	MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
	Vp_V_assert_sizes( v_lhs->dim(), rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1) );
	Vp_StMtV(v_lhs,alpha,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// v_lhs = op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_MtV(VectorWithOpMutable* v_lhs, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2)
{
	MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
	Vp_V_assert_sizes( v_lhs->dim(), rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1) );
	Vp_StMtV(v_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
// Level 3 BLAS

// M_lhs = alpha * op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_StMtM(MatrixWithOp* M_lhs, value_type alpha, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
						, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	assert_M_lhs(	  *M_lhs
					, rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
					, cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
	Mp_StMtM(M_lhs,alpha,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// M_lhs = op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_MtM(MatrixWithOp* M_lhs, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
						, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	assert_M_lhs(	  M_lhs
					, rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
					, cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
	Mp_StMtM(M_lhs,1.0,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0,0);
}

} // end namespace LinAlgOpPack


#endif // ABSTRACT_LIN_ALG_OP_PACK_DEF_H
