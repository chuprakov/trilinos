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

namespace LinAlgOpPack {

using BLAS_Cpp::rows;
using BLAS_Cpp::cols;

// Inject assert functions
using AbstractLinAlgPack::Vp_V_assert_compatibility;
using AbstractLinAlgPack::VopV_assert_compatibility;
using AbstractLinAlgPack::Mp_M_assert_compatibility;
using AbstractLinAlgPack::MopM_assert_compatibility;
using AbstractLinAlgPack::Vp_MtV_assert_compatibility;
using AbstractLinAlgPack::MtV_assert_compatibility;
using AbstractLinAlgPack::MtM_assert_compatibility;
using AbstractLinAlgPack::Mp_MtM_assert_compatibility;

// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Vectors

// v_lhs = V_rhs.
template <class V>
void assign(VectorWithOpMutable* v_lhs, const V& V_rhs) {
	Vp_V_assert_compatibility(v_lhs,V_rhs);
	(*v_lhs) = 0.0;
	Vp_V(v_lhs,V_rhs);
}

// v_lhs = alpha * V_rhs.
template <class V>
void V_StV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs) {
	Vp_V_assert_compatibility(v_lhs,V_rhs);
	(*v_lhs) = 0.0;
	Vp_StV(v_lhs,alpha,V_rhs);
}

// v_lhs = V1_rhs1 + V2_rhs2.
template <class V1, class V2>
void V_VpV(VectorWithOpMutable* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_compatibility(V1_rhs1,V2_rhs2);
	Vp_V_assert_compatibility(v_lhs,V1_rhs1);
	(*v_lhs) = 0.0;
	Vp_V(v_lhs,V1_rhs1);
	Vp_V(v_lhs,V2_rhs2);
}

// v_lhs = V_rhs1 - V_rhs2.
template <class V1, class V2>
void V_VmV(VectorWithOpMutable* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_compatibility(V1_rhs1,V2_rhs2);
	Vp_V_assert_compatibility(v_lhs,V1_rhs1);
	(*v_lhs) = 0.0;
	Vp_V(v_lhs,V1_rhs1);
	Vp_StV(v_lhs,-1.0,V2_rhs2);
}

// v_lhs = alpha * V_rhs1 + v_rhs2.
template <class V>
void V_StVpV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs1
	, const VectorWithOp& v_rhs2)
{
	VopV_assert_compatibility(V_rhs1,v_rhs2);
	(*v_lhs) = v_rhs2;
	Vp_StV(v_lhs,alpha,V_rhs1);
}

// ////////////////////////////////////////////////////////////////
// Level 2 BLAS

// v_lhs = alpha * op(M_rhs1) * V_rhs2.
template <class V>
void V_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2)
{
	Vp_MtV_assert_compatibility(v_lhs,M_rhs1,trans_rhs1,V_rhs2);
	Vp_StMtV(v_lhs,alpha,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// v_lhs = op(M_rhs1) * V_rhs2.
template <class V>
void V_MtV(VectorWithOpMutable* v_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2)
{
	Vp_MtV_assert_compatibility(v_lhs,M_rhs1,trans_rhs1,V_rhs2);
	Vp_StMtV(v_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

} // end namespace LinAlgOpPack

#endif // ABSTRACT_LIN_ALG_OP_PACK_DEF_H
