// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////////////
// TSFCoreVectorStdOps.hpp

#ifndef TSFCORE_VECTOR_STD_OPS_HPP
#define TSFCORE_VECTOR_STD_OPS_HPP

#include "TSFCoreVectorStdOpsDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
#include "RTOpPack_RTOpC.hpp"
#include "RTOp_ROp_norms.h"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_add_scalar.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_ele_wise_divide.h"
#include "RTOp_TOp_ele_wise_prod.h"
#include "RTOp_TOp_random_vector.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOp_TOp_sign.h"
#include "Teuchos_arrayArg.hpp"
#include "Teuchos_TestForException.hpp"

// RAB: 3/24/2004: ToDo: Replace all of these C RTOp subclases with
// fully templated direct subclasses of RTOpT<>!

// Reduction operations

template<class Scalar>
Scalar TSFCore::sum( const Vector<Scalar>& v_rhs )
{
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(sum_op,1,vecs,0,NULL,&*sum_targ);
	return sum_op(*sum_targ);
}

template<class Scalar>
Scalar TSFCore::norm_1( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC norm_1_op;
	TEST_FOR_EXCEPTION( 0!=RTOp_ROp_norm_1_construct(&norm_1_op.op()),std::logic_error,"Error!");
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(norm_1_op,1,vecs,0,NULL,&*norm_1_targ);
	return RTOp_ROp_norm_1_val(norm_1_op(*norm_1_targ));
}

template<class Scalar>
Scalar TSFCore::norm_2( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC norm_2_op;
	TEST_FOR_EXCEPTION( 0!=RTOp_ROp_norm_2_construct(&norm_2_op.op()),std::logic_error,"Error!");
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_2_targ = norm_2_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(norm_2_op,1,vecs,0,NULL,&*norm_2_targ);
	return RTOp_ROp_norm_2_val(norm_2_op(*norm_2_targ));
}

template<class Scalar>
Scalar TSFCore::norm_inf( const Vector<Scalar>& v_rhs )
{
	RTOpPack::RTOpC norm_inf_op;
	TEST_FOR_EXCEPTION( 0!=RTOp_ROp_norm_inf_construct(&norm_inf_op.op()),std::logic_error,"Error!");
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_inf_targ = norm_inf_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(norm_inf_op,1,vecs,0,NULL,&*norm_inf_targ);
	return RTOp_ROp_norm_inf_val(norm_inf_op(*norm_inf_targ));
}

template<class Scalar>
Scalar TSFCore::dot( const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2 )
{
	RTOpPack::RTOpC dot_prod_op;
	TEST_FOR_EXCEPTION(0!=RTOp_ROp_dot_prod_construct(&dot_prod_op.op()),std::logic_error,"Error!");
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> dot_prod_targ = dot_prod_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs1, &v_rhs2 };
	applyOp<Scalar>(dot_prod_op,2,vecs,0,NULL,&*dot_prod_targ);
	return RTOp_ROp_dot_prod_val(dot_prod_op(*dot_prod_targ));
}

template<class Scalar>
Scalar TSFCore::get_ele( const Vector<Scalar>& v, Index i )
{
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v };
	applyOp<Scalar>(sum_op,1,vecs,0,NULL,&*sum_targ,i,1);
	return sum_op(*sum_targ);
}

// Transformation operations

template<class Scalar>
void TSFCore::set_ele( Index i, Scalar alpha, Vector<Scalar>* v )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v==NULL,std::logic_error,"set_ele(...), Error!");
#endif
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
	Vector<Scalar>* targ_vecs[] = { v };
	applyOp<Scalar>(assign_scalar_op,0,NULL,1,targ_vecs,NULL,i,1);
}

template<class Scalar>
void TSFCore::assign( Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
	Vector<Scalar>* targ_vecs[] = { v_lhs };
	applyOp<Scalar>(assign_scalar_op,0,NULL,1,targ_vecs,NULL);
}

template<class Scalar>
void TSFCore::assign( Vector<Scalar>* v_lhs, const Vector<Scalar>& v_rhs )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::RTOpC assign_vectors_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_assign_vectors_construct(&assign_vectors_op.op())
    ,std::logic_error,"Error!" );
	const Vector<Scalar>* vecs[]      = { &v_rhs };
	Vector<Scalar>*       targ_vecs[] = { v_lhs  };
	applyOp<Scalar>(assign_vectors_op,1,vecs,1,targ_vecs,NULL);
}

template<class Scalar>
void TSFCore::Vp_S( Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
	RTOpPack::RTOpC add_scalar_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_add_scalar_construct(alpha,&add_scalar_op.op())
    ,std::logic_error,"Error!" );
	Vector<Scalar>* targ_vecs[] = { v_lhs };
	applyOp<Scalar>(add_scalar_op,0,NULL,1,targ_vecs,NULL);
}

template<class Scalar>
void TSFCore::Vt_S(
	Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
	if( alpha == 0.0 ) {
		assign(v_lhs,0.0);
	}
	else if( alpha != 1.0 ) {
		RTOpPack::RTOpC  scale_vector_op;
    TEST_FOR_EXCEPTION(
      0!=RTOp_TOp_scale_vector_construct(alpha,&scale_vector_op.op())
      ,std::logic_error,"Error!" );
		Vector<Scalar>* targ_vecs[] = { v_lhs };
		applyOp<Scalar>(scale_vector_op,0,NULL,1,targ_vecs,NULL);
	}
}

template<class Scalar>
void TSFCore::Vp_StV( Vector<Scalar>* v_lhs, const Scalar& alpha, const Vector<Scalar>& v_rhs )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vp_StV(...), Error!");
#endif
	RTOpPack::RTOpC axpy_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_axpy_construct(alpha,&axpy_op.op())
    ,std::logic_error,"Error!" );
	const Vector<Scalar>* vecs[]      = { &v_rhs };
	Vector<Scalar>*       targ_vecs[] = { v_lhs  };
	applyOp<Scalar>(axpy_op,1,vecs,1,targ_vecs,NULL);
}

template<class Scalar>
void TSFCore::ele_wise_prod(
	const Scalar& alpha, const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2
	,Vector<Scalar>* v_lhs
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_prod(...), Error");
#endif
	RTOpPack::RTOpC ele_wise_prod_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_ele_wise_prod_construct(alpha,&ele_wise_prod_op.op())
    ,std::logic_error,"Error!" );
	const Vector<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
	Vector<Scalar>*       targ_vecs[] = { v_lhs };
	applyOp<Scalar>(ele_wise_prod_op,2,vecs,1,targ_vecs,NULL);
}

template<class Scalar>
void TSFCore::ele_wise_divide(
	const Scalar& alpha, const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2
	,Vector<Scalar>* v_lhs
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_divide(...), Error");
#endif
	RTOpPack::RTOpC ele_wise_divide_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_ele_wise_divide_construct(alpha,&ele_wise_divide_op.op())
    ,std::logic_error,"Error!" );
	const Vector<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
	Vector<Scalar>*       targ_vecs[] = { v_lhs };
	applyOp<Scalar>(ele_wise_divide_op,2,vecs,1,targ_vecs,NULL);
}

template<class Scalar>
void TSFCore::seed_randomize( unsigned int s )
{
	srand(s);
}

template<class Scalar>
void TSFCore::randomize( Scalar l, Scalar u, Vector<Scalar>* v )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v==NULL,std::logic_error,"Vt_S(...), Error");
#endif
	RTOpPack::RTOpC random_vector_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_random_vector_construct(l,u,&random_vector_op.op())
    ,std::logic_error,"Error!" );
	Vector<Scalar>* targ_vecs[] = { v };
	applyOp<Scalar>(random_vector_op,0,NULL,1,targ_vecs,NULL);
}

#endif // TSFCORE_VECTOR_STD_OPS_HPP
