// //////////////////////////////////////////////////////////////////////////////////
// LinAlgOpPackDecl.h
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

#ifndef ABSTRACT_LIN_ALG_OP_PACK_DECL_H
#define ABSTRACT_LIN_ALG_OP_PACK_DECL_H

#include "AbstractLinAlgPackTypes.h"
#include "MatrixWithOpNonsingular.h"
#include "VectorWithOpMutable.h"
#include "VectorStdOps.h"
#include "AbstractLinAlgPackAssertOp.h"

namespace LinAlgOpPack {

///
typedef AbstractLinAlgPack::size_type  size_type;
///
typedef AbstractLinAlgPack::value_type value_type;

///
using AbstractLinAlgPack::VectorSpace;
///
using AbstractLinAlgPack::VectorWithOp;
///
using AbstractLinAlgPack::VectorWithOpMutable;
///
using AbstractLinAlgPack::MatrixWithOp;
///
using AbstractLinAlgPack::MatrixNonsingular;
///
using AbstractLinAlgPack::MatrixWithOpNonsingular;

// Inject names of base linear algebra functions for LinAlgPack.
// Note that this is neccesary in MS VC++ 5.0 because
// it does not perform name lookups properly but it
// is not adverse to the standard so it is a portable
// fix.
///
using AbstractLinAlgPack::sum;
///
using AbstractLinAlgPack::dot;
///
using AbstractLinAlgPack::Vp_S;
///
using AbstractLinAlgPack::Vt_S;
///
using AbstractLinAlgPack::Vp_StV;
///
using AbstractLinAlgPack::Vp_StMtV;
///
using AbstractLinAlgPack::Mt_S;
///
using AbstractLinAlgPack::Mp_StM;
///
using AbstractLinAlgPack::Mp_StMtM;
///
using AbstractLinAlgPack::V_InvMtV;
///
using AbstractLinAlgPack::M_StInvMtM;

/** \defgroup LinAlgOpPack_grp Default linear algebra implementation operations.
  *
  * These are template functions that can be used to perform simpler
  * linear algebra operations given more elaborate ones.  The idea is that
  * for each combination of vector and matrix types, the BLAS like operations
  * must be provided and then these template functions provide related
  * linear algebra operations.  The user can override these default implementations
  * by defining the exact functions himself.
  *
  * Warning!
  * In general it is not allowed for the lhs argument to be used in the rhs expression.
  * Concidering aliasing would make the operations much more complicated to implement.
  * So unless you are sure that it is okay, do not use a vector or matrix object
  * in both the lhs and rhs expressions.
  */
//@{

/** @name Level 1 BLAS for Vectors
  *
  * For these functions to work for the type V the following function must
  * be defined:
  \code
  // v_lhs += alpha * V_rhs
  void Vp_StV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs);
  \endcode
  * The rest of these level 1 BLAS functions implement the variations.
  */
//@{

///
/** v_lhs += V_rhs.
  *
  * Calls: <tt>Vp_StV(v_lhs,1.0,V_rhs);</tt>
  */
template <class V>
void Vp_V(VectorWithOpMutable* v_lhs, const V& V_rhs);

///
/** v_lhs = V_rhs.
  *
  * Calls: <tt>Vp_V(v_lhs,V_rhs);</tt>
  */
template <class V>
void assign(VectorWithOpMutable* v_lhs, const V& V_rhs);

///
/** v_lhs = alpha * V_rhs.
  *
  * Calls: <tt>Vp_StV(v_lhs,alpha,V_rhs);</tt>
  */
template <class V>
void V_StV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs);

///
/** v_lhs = - V_rhs.
  *
  * Calls: <tt>V_StV(v_lhs,-1.0,V_rhs);</tt>
  */
template <class V>
void V_mV(VectorWithOpMutable* v_lhs, const V& V_rhs);

/// 
/** v_lhs = V1_rhs1 + V2_rhs2.
  *
  * Calls: <tt>Vp_V(v_lhs,V1_rhs1); Vp_V(v_lhs,V1_rhs2);</tt>
  */
template <class V1, class V2>
void V_VpV(VectorWithOpMutable* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2);

///
/** v_lhs = V_rhs1 - V_rhs2.
  *
  * Calls: <tt>Vp_V(v_lhs,V1_rhs1); Vp_StV(v_lhs,-1.0,V2_rhs2);</tt>
  */
template <class V1, class V2>
void V_VmV(VectorWithOpMutable* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2);

///
/** v_lhs = alpha * V_rhs1 + vs_rhs2.
  *
  * Calls: <tt>Vp_StV(v_lhs,alpha,V_rhs1);</tt>
  */
template <class V>
void V_StVpV(VectorWithOpMutable* v_lhs, value_type alpha, const V& V_rhs1
	, const VectorWithOp& vs_rhs2);

//		end Level 1 BLAS for Vectors
//@}

// //////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
/** @name Level 1 BLAS for Matrices
  *
  * For these functions to work for the type M the following function must
  * be defined:
  *
  * // M_lhs += alpha * op(M_rhs)	\\
  * void Mp_StM(MatrixWithOp* v_lhs, value_type alpha, const V& V_rhs, BLAS_Cpp::Transp);
  *
  * The rest of these level 1 BLAS functions implement the variations.
  */
//@{

///
/** M_lhs += op(M_rhs).
  *
  * Calls: <tt>Mp_StM(M_lhs,1.0,M_rhs,trans_rhs);</tt>
  */
void Mp_M(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs);

//		end += operations
//@}

///
/** M_lhs = op(M_rhs).
  *
  * Calls: <tt>Mp_M(M_lhs,M_rhs,trans_rhs);</tt>
  */
void assign(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs);

///
/** gm_lhs = alpha * M_rhs.
  *
  * Calls: <tt>Mp_StM(&(*gm_lhs)(),alpha,M_rhs,trans_rhs);</tt>
  */
void M_StM(MatrixWithOp* M_lhs, value_type alpha, const MatrixWithOp& M_rhs
	, BLAS_Cpp::Transp trans_rhs);

///
/** gm_lhs = - op(M_rhs).
  *
  * Calls: <tt>M_StM(&(*gm_lhs)(),-1.0,M_rhs,trans_rhs);</tt>
  */
void M_mM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs) ;

/// 
/** M_lhs = op(M_rhs1) + op(M_rhs2).
  *
  * Calls: <tt>Mp_M(M_lhs,M_rhs1,trans_rhs1); Mp_M(M_lhs,M_rhs2,trans_rhs2);</tt>
  */
void M_MpM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2);

///
/** M_lhs = op(M_rhs1) - op(M_rhs2).
  *
  * Calls: <tt>Mp_M(M_lhs,M_rhs1,trans_rhs1); Mp_StM(M_lhs,-1.0,M_rhs2,trans_rhs2);</tt>
  */
void M_MmM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2);

///
/** M_lhs = alpha * op(M_rhs1) + op(gms_rhs2).
  *
  * Calls: <tt>Mp_StM(M_lhs,alpha,M_rhs1,trans_rhs1);</tt>
  */
void M_StMpM(MatrixWithOp* M_lhs, value_type alpha, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& gms_rhs2, BLAS_Cpp::Transp trans_rhs2);

//		end Level 1 BLAS for Matrices
//@}

// //////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////// 
/** @name Level 2 BLAS
  *
  * These operations implement variations on the Level-2 BLAS operation:\\
  *
  * v_lhs = alpha * op(M_rhs1) * V_rhs2 + beta * v_lhs
  */
//@{

///
/** v_lhs += op(M_rhs1) * V_rhs2.
  *
  * Calls: <tt>Vp_StMtV(v_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2);</tt>
  */
template <class V>
void Vp_MtV(VectorWithOpMutable* v_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2);

///
/** v_lhs = alpha * op(M_rhs1) * V_rhs2.
  *
  * Calls: <tt>Vp_StMtV(v_lhs,alpha,M_rhs1,trans_rhs1,V_rhs2);</tt>
  */
template <class V>
void V_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2);

///
/** v_lhs = op(M_rhs1) * V_rhs2.
  *
  * Calls: <tt>Vp_MtV(v_lhs,M_rhs1,trans_rhs1,V_rhs2);</tt>
  */
template <class V>
void V_MtV(VectorWithOpMutable* v_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2);

//		end Level 2 BLAS
//@}

// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
/** @name Level 3 BLAS
  *
  * These operations are based on the Level-3 BLAS operation:
  *
  * M_lhs = alpha * op(M_rhs1) * op(M_rhs2) + beta * M_lhs
  */
//@{

///
/** M_lhs += op(M_rhs1) * op(M_rhs2).
  *
  * Calls: <tt>Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);</tt>
  */
void Mp_MtM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2);

///
/** M_lhs = op(M_rhs1) * op(M_rhs2) + beta * gms_rhs.
  *
  * Calls: <tt>Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,beta);</tt>
  */
void Mp_MtM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2, value_type beta);

///
/** M_lhs = alpha * op(M_rhs1) * op(M_rhs2).
  *
  * Calls: <tt>Mp_StMtM(M_lhs,alpha,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);</tt>
  */
void M_StMtM(MatrixWithOp* M_lhs, value_type alpha, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2);

///
/** M_lhs = op(M_rhs1) * op(M_rhs2).
  *
  * Calls: <tt>Mp_MtM(M_lhs,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);</tt>
  */
void M_MtM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2);

//		end Level 3 BLAS
//@}

//		end Default Linear Algebra Implementations
//@}

// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////
// Inline definitions

// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Vectors

// v_lhs += V_rhs.
template <class V>
inline
void Vp_V(VectorWithOpMutable* v_lhs, const V& V_rhs) {
	Vp_StV(v_lhs,1.0,V_rhs);
}

// v_lhs = - V_rhs.
template <class V>
inline
void V_mV(VectorWithOpMutable* v_lhs, const V& V_rhs) {
	V_StV(v_lhs,-1.0,V_rhs);
}

// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Matrices

// M_lhs += op(M_rhs).
inline
void Mp_M(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	Mp_StM(M_lhs,1.0,M_rhs,trans_rhs);
}

// M_lhs = - op(M_rhs).
inline
void M_mM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs) {
	M_StM(M_lhs,-1.0,M_rhs,trans_rhs);
}

// /////////////////////////////////////////////////////////////////////// 
// Level 2 BLAS

// v_lhs += op(M_rhs1) * V_rhs2.
template <class V>
inline
void Vp_MtV(VectorWithOpMutable* v_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2)
{
	Vp_StMtV(v_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2);
}

// v_lhs = op(M_rhs1) * V_rhs2 + beta * v_lhs.
template <class V>
inline
void Vp_MtV(VectorWithOpMutable* v_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const V& V_rhs2, value_type beta)
{
	Vp_StMtV(v_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2,beta);
}

// //////////////////////////////////////////////////////////////////////////////
// Level 3 BLAS

// M_lhs += op(M_rhs1) * op(M_rhs2).
inline
void Mp_MtM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
}

// M_lhs = op(M_rhs1) * op(M_rhs2) + beta * gms_rhs.
inline
void Mp_MtM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,beta);
}

// M_lhs = inv(op(M_rhs1)) * op(M_rhs2)
inline
void M_InvMtM(MatrixWithOp* M_lhs, const MatrixNonsingular& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{
	M_StInvMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
}

} // end namespace LinAlgOpPack

#endif	// ABSTRACT_LIN_ALG_OP_PACK_DECL_H
