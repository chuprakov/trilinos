// ////////////////////////////////////////////////////////////////////
// VectorStdOps.hpp

#ifndef TSFCORE_VECTOR_STD_OPS_HPP
#define TSFCORE_VECTOR_STD_OPS_HPP

#include <assert.h>

#include "TSFCoreVectorStdOpsDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "RTOp_ROp_sum.h"
#include "RTOp_ROp_norms.h"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_add_scalar.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_ele_wise_divide.h"
#include "RTOp_TOp_ele_wise_prod.h"
#include "RTOp_TOp_random_vector.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOp_TOp_sign.h"
#include "RTOpCppC.hpp"
#include "ThrowException.hpp"

// Reduction operations

template<class Scalar>
Scalar TSFCore::sum( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC          sum_op;
	RTOpPack::ReductTarget   sum_targ;
	if(0>RTOp_ROp_sum_construct(&sum_op.op() )) assert(0);
	sum_op.reduct_obj_create(&sum_targ);
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp(sum_op,1,vecs,0,NULL,sum_targ.obj() );
	return RTOp_ROp_sum_val(sum_targ.obj());
}

template<class Scalar>
Scalar TSFCore::norm_1( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC          norm_1_op;
	RTOpPack::ReductTarget   norm_1_targ;
	if(0>RTOp_ROp_norm_1_construct(&norm_1_op.op() )) assert(0);
	norm_1_op.reduct_obj_create(&norm_1_targ);
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp(norm_1_op,1,vecs,0,NULL,norm_1_targ.obj() );
	return RTOp_ROp_norm_1_val(norm_1_targ.obj());
}

template<class Scalar>
Scalar TSFCore::norm_2( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC          norm_2_op;
	RTOpPack::ReductTarget   norm_2_targ;
	if(0>RTOp_ROp_norm_2_construct(&norm_2_op.op() )) assert(0);
	norm_2_op.reduct_obj_create(&norm_2_targ);
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp(norm_2_op,1,vecs,0,NULL,norm_2_targ.obj() );
	return RTOp_ROp_norm_2_val(norm_2_targ.obj());
}

template<class Scalar>
Scalar TSFCore::norm_inf( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC          norm_inf_op;
	RTOpPack::ReductTarget   norm_inf_targ;
	if(0>RTOp_ROp_norm_inf_construct(&norm_inf_op.op() )) assert(0);
	norm_inf_op.reduct_obj_create(&norm_inf_targ);
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp(norm_inf_op,1,vecs,0,NULL,norm_inf_targ.obj() );
	return RTOp_ROp_norm_inf_val(norm_inf_targ.obj());
}

template<class Scalar>
Scalar TSFCore::dot( const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2 )
{
	RTOpPack::RTOpC          dot_prod_op;
	RTOpPack::ReductTarget   dot_prod_targ;
	if(0>RTOp_ROp_dot_prod_construct(&dot_prod_op.op())) assert(0);
	dot_prod_op.reduct_obj_create(&dot_prod_targ);
	const Vector<Scalar>* vecs[] = { &v_rhs1, &v_rhs2 };
	applyOp<Scalar>(dot_prod_op,2,vecs,0,NULL,dot_prod_targ.obj() );
	return RTOp_ROp_dot_prod_val(dot_prod_targ.obj());
}

template<class Scalar>
Scalar TSFCore::get_ele( const Vector<Scalar>& v, Index i )
{
	RTOpPack::RTOpC          sum_op;
	RTOpPack::ReductTarget   sum_targ;
	if(0>RTOp_ROp_sum_construct(&sum_op.op() )) assert(0);
	sum_op.reduct_obj_create(&sum_targ);
	const Vector<Scalar>* vecs[] = { &v };
	applyOp<Scalar>(sum_op,1,vecs,0,NULL,sum_targ.obj(),i,1);
	return RTOp_ROp_sum_val(sum_targ.obj());
}

// Transformation operations

template<class Scalar>
void TSFCore::set_ele( Index i, Scalar alpha, Vector<Scalar>* v )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v==NULL,std::logic_error,"set_ele(...), Error!");
#endif
	RTOpPack::RTOpC  assign_scalar_op;
	if(0>RTOp_TOp_assign_scalar_construct(alpha,&assign_scalar_op.op())) assert(0);
	Vector<Scalar>* targ_vecs[] = { v };
	applyOp<Scalar>(assign_scalar_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL,i,1);
}

template<class Scalar>
void TSFCore::assign( Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::RTOpC  assign_scalar_op;
	if(0>RTOp_TOp_assign_scalar_construct(alpha,&assign_scalar_op.op())) assert(0);
	Vector<Scalar>* targ_vecs[] = { v_lhs };
	applyOp<Scalar>(assign_scalar_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::assign( Vector<Scalar>* v_lhs, const Vector<Scalar>& v_rhs )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::RTOpC assign_vectors_op;
	if(0>RTOp_TOp_assign_vectors_construct(&assign_vectors_op.op())) assert(0);
	const Vector<Scalar>* vecs[]      = { &v_rhs };
	Vector<Scalar>*       targ_vecs[] = { v_lhs  };
	applyOp<Scalar>(assign_vectors_op,1,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::Vp_S( Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
	RTOpPack::RTOpC  add_scalar_op;
	if(0>RTOp_TOp_add_scalar_construct(alpha,&add_scalar_op.op())) assert(0);
	Vector<Scalar>* targ_vecs[] = { v_lhs };
	applyOp<Scalar>(add_scalar_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::Vt_S(
	Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
	if( alpha == 0.0 ) {
		assign(v_lhs,0.0);
	}
	else if( alpha != 1.0 ) {
		RTOpPack::RTOpC  scale_vector_op;
		if(0>RTOp_TOp_scale_vector_construct(alpha,&scale_vector_op.op())) assert(0);
		Vector<Scalar>* targ_vecs[] = { v_lhs };
		applyOp<Scalar>(scale_vector_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
	}
}

template<class Scalar>
void TSFCore::Vp_StV( Vector<Scalar>* v_lhs, const Scalar& alpha, const Vector<Scalar>& v_rhs )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"Vp_StV(...), Error!");
#endif
	RTOpPack::RTOpC axpy_op;
	if(0>RTOp_TOp_axpy_construct(alpha,&axpy_op.op())) assert(0);
	const Vector<Scalar>* vecs[]      = { &v_rhs };
	Vector<Scalar>*       targ_vecs[] = { v_lhs  };
	applyOp<Scalar>(axpy_op,1,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::ele_wise_prod(
	const Scalar& alpha, const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2
	,Vector<Scalar>* v_lhs
	)
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_prod(...), Error");
#endif
	RTOpPack::RTOpC ele_wise_prod_op;
	if(0>RTOp_TOp_ele_wise_prod_construct(alpha,&ele_wise_prod_op.op())) assert(0);
	const Vector<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
	Vector<Scalar>*       targ_vecs[] = { v_lhs };
	applyOp<Scalar>(ele_wise_prod_op,2,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::ele_wise_divide(
	const Scalar& alpha, const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2
	,Vector<Scalar>* v_lhs
	)
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_divide(...), Error");
#endif
	RTOpPack::RTOpC ele_wise_divide_op;
	if(0>RTOp_TOp_ele_wise_divide_construct(alpha,&ele_wise_divide_op.op())) assert(0);
	const Vector<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
	Vector<Scalar>*       targ_vecs[] = { v_lhs };
	applyOp<Scalar>(ele_wise_divide_op,2,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::seed_random_vector_generator( unsigned int s )
{
	srand(s);
}

template<class Scalar>
void TSFCore::random_vector( Scalar l, Scalar u, Vector<Scalar>* v )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v==NULL,std::logic_error,"Vt_S(...), Error");
#endif
	RTOpPack::RTOpC  random_vector_op;
	if(0>RTOp_TOp_random_vector_construct(l,u,&random_vector_op.op())) assert(0);
	Vector<Scalar>* targ_vecs[] = { v };
	applyOp<Scalar>(random_vector_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

#endif // TSFCORE_VECTOR_STD_OPS_HPP
