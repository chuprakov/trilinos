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

#ifndef LIN_ALG_OP_PACK_DEF_H
#define LIN_ALG_OP_PACK_DEF_H

#include "LinAlgOpPackDecl.h"	// also includes some inline function definitions
#include "LinAlgPackAssertOp.h"
#include "GenMatrixClass.h"

namespace LinAlgOpPack {

using BLAS_Cpp::rows;
using BLAS_Cpp::cols;

// Inject assert functions
using LinAlgPack::assert_gms_lhs;
using LinAlgPack::Vp_V_assert_sizes;
using LinAlgPack::VopV_assert_sizes;
using LinAlgPack::Mp_M_assert_sizes;
using LinAlgPack::MopM_assert_sizes;
using LinAlgPack::Vp_MtV_assert_sizes;
using LinAlgPack::MtV_assert_sizes;
using LinAlgPack::MtM_assert_sizes;

// Inject names of base linear algebra functions for LinAlgPack.
// Note that this is neccesary in MS VC++ 5.0 because
// it does not perform name lookups properly but it
// is not adverse to the standard so it is a portable
// fix.
using LinAlgPack::assign;
using LinAlgPack::Vt_S;
using LinAlgPack::Vp_StV;
using LinAlgPack::Vp_StMtV;
using LinAlgPack::Mt_S;
using LinAlgPack::Mp_StM;
using LinAlgPack::Mp_StMtM;

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Vectors

// //////////////////////////////////////////////////////////////////////////////
// += operations 

// //////////////////////////////////////////////////////////////////////////////
// operations with Vector as lhs

// v_lhs = V_rhs.
template <class V>
void assign(Vector* v_lhs, const V& V_rhs) {
	v_lhs->resize(V_rhs.dim());
	(*v_lhs) = 0.0;
	Vp_V(&(*v_lhs)(),V_rhs);
}

// v_lhs = alpha * V_rhs.
template <class V>
void V_StV(Vector* v_lhs, value_type alpha, const V& V_rhs) {
	v_lhs->resize(V_rhs.dim());
	(*v_lhs) = 0.0;
	Vp_StV(&(*v_lhs)(),alpha,V_rhs);
}

// v_lhs = V1_rhs1 + V2_rhs2.
template <class V1, class V2>
void V_VpV(Vector* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
	v_lhs->resize(V1_rhs1.dim());
	(*v_lhs) = 0.0;
	VectorSlice vs_lhs(*v_lhs);
	Vp_V(&vs_lhs,V1_rhs1);
	Vp_V(&vs_lhs,V2_rhs2);
}


// v_lhs = V_rhs1 - V_rhs2.
template <class V1, class V2>
void V_VmV(Vector* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
	v_lhs->resize(V1_rhs1.dim());
	(*v_lhs) = 0.0;
	VectorSlice vs_lhs(*v_lhs);
	Vp_V(&vs_lhs,V1_rhs1);
	Vp_StV(&vs_lhs,-1.0,V2_rhs2);
}


// v_lhs = alpha * V_rhs1 + vs_rhs2.
template <class V>
void V_StVpV(Vector* v_lhs, value_type alpha, const V& V_rhs1
	, const VectorSlice& vs_rhs2)
{
	VopV_assert_sizes(V_rhs1.dim(),vs_rhs2.dim());
	(*v_lhs) = vs_rhs2;
	Vp_StV(&(*v_lhs)(),alpha,V_rhs1);
}

// ///////////////////////////////////////////////////////////////////////////
// operations with VectorSlice as lhs

// vs_lhs = V_rhs.
template <class V>
void assign(VectorSlice* vs_lhs, const V& V_rhs) {
	Vp_V_assert_sizes( vs_lhs->dim(), V_rhs.dim() );
	(*vs_lhs) = 0.0;
	Vp_V(vs_lhs,V_rhs);
}

// vs_lhs = alpha * V_rhs.
template <class V>
void V_StV(VectorSlice* vs_lhs, value_type alpha, const V& V_rhs) {
	Vp_V_assert_sizes( vs_lhs->dim(), V_rhs.dim() );
	(*vs_lhs) = 0.0;
	Vp_StV(vs_lhs,alpha,V_rhs);
}

// vs_lhs = V1_rhs1 + V2_rhs2.
template <class V1, class V2>
void V_VpV(VectorSlice* vs_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
	Vp_V_assert_sizes( vs_lhs->dim(), V1_rhs1.dim() );
	(*vs_lhs) = 0.0;
	Vp_V(vs_lhs,V1_rhs1);
	Vp_V(vs_lhs,V2_rhs2);
}

// vs_lhs = V_rhs1 - V_rhs2.
template <class V1, class V2>
void V_VmV(VectorSlice* vs_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
	VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
	Vp_V_assert_sizes( vs_lhs->dim(), V1_rhs1.dim() );
	(*vs_lhs) = 0.0;
	Vp_V(vs_lhs,V1_rhs1);
	Vp_StV(vs_lhs,-1.0,V2_rhs2);
}

// vs_lhs = alpha * V_rhs1 + vs_rhs2.
template <class V>
void V_StVpV(VectorSlice* vs_lhs, value_type alpha, const V& V_rhs1
	, const VectorSlice& vs_rhs2)
{
	VopV_assert_sizes(V_rhs1.dim(),vs_rhs2.dim());
	(*vs_lhs) = vs_rhs2;
	Vp_StV(vs_lhs,alpha,V_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Matrices

// //////////////////////////////////////////////////////////////////////////////
// += operations 


// //////////////////////////////////////////////////////////////////////////////
// operations with GenMatrix as lhs

// gm_lhs = op(M_rhs).
template <class M>
void assign(GenMatrix* gm_lhs, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	gm_lhs->resize(	 rows(M_rhs.rows(),M_rhs.cols(),trans_rhs)
					,cols(M_rhs.rows(),M_rhs.cols(),trans_rhs) );
	(*gm_lhs) = 0.0;
	Mp_StM(&(*gm_lhs)(),1.0,M_rhs,trans_rhs);
}

// gm_lhs = alpha * op(M_rhs).
template <class M>
void M_StM(GenMatrix* gm_lhs, value_type alpha, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	gm_lhs->resize(	 rows(M_rhs.rows(),M_rhs.cols(),trans_rhs)
					,cols(M_rhs.rows(),M_rhs.cols(),trans_rhs) );
	(*gm_lhs) = 0.0;
	Mp_StM(&(*gm_lhs)(),alpha,M_rhs,trans_rhs);
}

// gm_lhs = op(M1_rhs1) + op(M2_rhs2).
template <class M1, class M2>
void M_MpM(GenMatrix* gm_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
						,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
	gm_lhs->resize(	 rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
					,cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs2) );
	(*gm_lhs) = 0.0;
	GenMatrixSlice gms_lhs(*gm_lhs);
	Mp_M(&gms_lhs,M1_rhs1,trans_rhs1);
	Mp_M(&gms_lhs,M2_rhs2,trans_rhs2);
}

// gm_lhs = op(M_rhs1) - op(M_rhs2).
template <class M1, class M2>
void M_MmM(GenMatrix* gm_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
						,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
	gm_lhs->resize(	 rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
					,cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
	(*gm_lhs) = 0.0;
	GenMatrixSlice gms_lhs(*gm_lhs);
	Mp_M(&gms_lhs,M1_rhs1,trans_rhs1);
	Mp_StM(&gms_lhs,-1.0,M2_rhs2,trans_rhs2);
}

// gm_lhs = alpha * op(M_rhs1) + op(gms_rhs2).
template <class M>
void M_StMpM(GenMatrix* gm_lhs, value_type alpha, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const GenMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M_rhs1.rows(),M_rhs1.cols(),trans_rhs1
						,gms_rhs2.rows(),gms_rhs2.cols(),trans_rhs2);
	assign(gm_lhs,gms_rhs2,trans_rhs2);
	Mp_StM(&(*gm_lhs)(),alpha,M_rhs1,trans_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// operations with GenMatrixSlice as lhs

// gms_lhs = op(M_rhs).
template <class M>
void assign(GenMatrixSlice* gms_lhs, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, M_rhs.rows(), M_rhs.cols(), trans_rhs	);
	(*gms_lhs) = 0.0;
	Mp_StM(gms_lhs,1.0,M_rhs,trans_rhs);
}

// gms_lhs = alpha * op(M_rhs).
template <class M>
void M_StM(GenMatrixSlice* gms_lhs, value_type alpha, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, M_rhs.rows(), M_rhs.cols(), trans_rhs	);
	(*gms_lhs) = 0.0;
	Mp_StM(gms_lhs,alpha,M_rhs,trans_rhs);
}

// gms_lhs = op(M1_rhs1) + op(M2_rhs2).
template <class M1, class M2>
void M_MpM(GenMatrixSlice* gms_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
						,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
	assert_gms_lhs(*gms_lhs, rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
						   , cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
	(*gms_lhs) = 0.0;
	Mp_M(gms_lhs,M1_rhs1,trans_rhs1);
	Mp_M(gms_lhs,M2_rhs2,trans_rhs2);
}

// gms_lhs = op(M_rhs1) - op(M_rhs2).
template <class M1, class M2>
void M_MmM(GenMatrixSlice* gms_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
						,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
	assert_gms_lhs(*gms_lhs, rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
						   , cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
	(*gms_lhs) = 0.0;
	Mp_M(gms_lhs,M1_rhs1,trans_rhs1);
	Mp_StM(gms_lhs,-1.0,M2_rhs2,trans_rhs2);
}

// gms_lhs = alpha * op(M_rhs1) + op(gms_rhs2).
template <class M>
void M_StMpM(GenMatrixSlice* gms_lhs, value_type alpha, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const GenMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MopM_assert_sizes(	 M_rhs1.rows(),M_rhs1.cols(),trans_rhs1
						,gms_rhs2.rows(),gms_rhs2.cols(),trans_rhs2);
	assign(gms_lhs,gms_rhs2,trans_rhs2);
	Mp_StM(gms_lhs,alpha,M_rhs1,trans_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////// /////
// Level 2 BLAS

// //////////////////////////////////////////////////////////////////////////////
// += operations

// //////////////////////////////////////////////////////////////////////////////
// operations with Vector as lhs

// v_lhs = alpha * op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_StMtV(Vector* v_lhs, value_type alpha, const M& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2)
{
	MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
	v_lhs->resize(rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1));
	Vp_StMtV(&(*v_lhs)(),alpha,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// v_lhs = op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_MtV(Vector* v_lhs, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2)
{
	MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
	v_lhs->resize(rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1));
	Vp_StMtV(&(*v_lhs)(),1.0,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// operations with VectorSlice as lhs

// vs_lhs = alpha * op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_StMtV(VectorSlice* vs_lhs, value_type alpha, const M& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2)
{
	MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
	Vp_V_assert_sizes( vs_lhs->dim(), rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1) );
	Vp_StMtV(vs_lhs,alpha,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// vs_lhs = op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_MtV(VectorSlice* vs_lhs, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2)
{
	MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
	Vp_V_assert_sizes( vs_lhs->dim(), rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1) );
	Vp_StMtV(vs_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
// Level 3 BLAS

// //////////////////////////////////////////////////////////////////////////////
// += operations 

// //////////////////////////////////////////////////////////////////////////////
// = operations with GenMatrix as lhs

// gm_lhs = alpha * op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_StMtM(GenMatrix* gm_lhs, value_type alpha, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
						, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	gm_lhs->resize(	  rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
					, cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
	Mp_StMtM(&(*gm_lhs)(),alpha,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// gm_lhs = op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_MtM(GenMatrix* gm_lhs, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
						, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	gm_lhs->resize(	  rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
					, cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
	Mp_StMtM(&(*gm_lhs)(),1.0,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// = operations with GenMatrixSlice as lhs

// gms_lhs = alpha * op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
						, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	assert_gms_lhs(	  *gms_lhs
					, rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
					, cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
	Mp_StMtM(gms_lhs,alpha,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// gms_lhs = op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_MtM(GenMatrixSlice* gms_lhs, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
						, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	assert_gms_lhs(	  gms_lhs
					, rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
					, cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
	Mp_StMtM(gms_lhs,1.0,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0,0);
}

} // end namespace LinAlgOpPack


#endif // LIN_ALG_OP_PACK_DEF_H
