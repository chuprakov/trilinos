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
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNorm2.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_ROpWeightedNorm2.hpp"
#include "RTOpPack_TOpAbs.hpp"
#include "RTOpPack_TOpAddScalar.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
#include "RTOpPack_TOpAssignVectors.hpp"
#include "RTOpPack_TOpAXPY.hpp"
#include "RTOpPack_TOpEleWiseDivide.hpp"
#include "RTOpPack_TOpEleWiseProd.hpp"
#include "RTOpPack_TOpLinearCombination.hpp"
#include "RTOpPack_TOpScaleVector.hpp"
#include "RTOpPack_TOpReciprocal.hpp"
#include "RTOpPack_TOpRandomize.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_arrayArg.hpp"

// Reduction operations

template<class Scalar>
Scalar TSFCore::sum( const Vector<Scalar>& v_rhs )
{
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(sum_op,1,vecs,0,(Vector<Scalar>**)NULL,&*sum_targ);
	return sum_op(*sum_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TSFCore::norm_1( const Vector<Scalar>& v_rhs )
{
  RTOpPack::ROpNorm1<Scalar> norm_1_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(norm_1_op,1,vecs,0,(Vector<Scalar>**)NULL,&*norm_1_targ);
	return norm_1_op(*norm_1_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TSFCore::norm_2( const Vector<Scalar>& v_rhs )
{
  RTOpPack::ROpNorm2<Scalar> norm_2_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_2_targ = norm_2_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(norm_2_op,1,vecs,0,(Vector<Scalar>**)NULL,&*norm_2_targ);
	return norm_2_op(*norm_2_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TSFCore::norm_2( const Vector<Scalar>& w, const Vector<Scalar>& v )
{
  RTOpPack::ROpWeightedNorm2<Scalar> wght_norm_2_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> wght_norm_2_targ = wght_norm_2_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &w, &v };
	applyOp<Scalar>(wght_norm_2_op,2,vecs,0,(Vector<Scalar>**)NULL,&*wght_norm_2_targ);
	return wght_norm_2_op(*wght_norm_2_targ);
}

template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TSFCore::norm_inf( const Vector<Scalar>& v_rhs )
{
  RTOpPack::ROpNormInf<Scalar> norm_inf_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> norm_inf_targ = norm_inf_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs };
	applyOp<Scalar>(norm_inf_op,1,vecs,0,(Vector<Scalar>**)NULL,&*norm_inf_targ);
	return norm_inf_op(*norm_inf_targ);
}

template<class Scalar>
Scalar TSFCore::dot( const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2 )
{
  RTOpPack::ROpDotProd<Scalar> dot_prod_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> dot_prod_targ = dot_prod_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v_rhs1, &v_rhs2 };
	applyOp<Scalar>(dot_prod_op,2,vecs,0,(Vector<Scalar>**)NULL,&*dot_prod_targ);
  return dot_prod_op(*dot_prod_targ);
}

template<class Scalar>
Scalar TSFCore::get_ele( const Vector<Scalar>& v, Index i )
{
  RTOpPack::ROpSum<Scalar> sum_op;
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> sum_targ = sum_op.reduct_obj_create();
	const Vector<Scalar>* vecs[] = { &v };
	applyOp<Scalar>(sum_op,1,vecs,0,(Vector<Scalar>**)NULL,&*sum_targ,i,1,0);
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
	applyOp<Scalar>(assign_scalar_op,0,(const Vector<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL,i,1,0);
}

template<class Scalar>
void TSFCore::assign( Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
	Vector<Scalar>* targ_vecs[] = { v_lhs };
	applyOp<Scalar>(assign_scalar_op,0,(const Vector<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::assign( Vector<Scalar>* v_lhs, const Vector<Scalar>& v_rhs )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAssignVectors<Scalar> assign_vectors_op;
	const Vector<Scalar>* vecs[]      = { &v_rhs };
	Vector<Scalar>*       targ_vecs[] = { v_lhs  };
	applyOp<Scalar>(assign_vectors_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::Vp_S( Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
  RTOpPack::TOpAddScalar<Scalar> add_scalar_op(alpha);
	Vector<Scalar>* targ_vecs[] = { v_lhs };
	applyOp<Scalar>(add_scalar_op,0,(const Vector<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::Vt_S(
	Vector<Scalar>* v_lhs, const Scalar& alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
	if( alpha == Teuchos::ScalarTraits<Scalar>::zero() ) {
		assign(v_lhs,Teuchos::ScalarTraits<Scalar>::zero());
	}
	else if( alpha != Teuchos::ScalarTraits<Scalar>::one() ) {
    RTOpPack::TOpScaleVector<Scalar> scale_vector_op(alpha);
		Vector<Scalar>* targ_vecs[] = { v_lhs };
		applyOp<Scalar>(scale_vector_op,0,(const Vector<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
	}
}

template<class Scalar>
void TSFCore::V_StV( Vector<Scalar>* y, const Scalar& alpha, const Vector<Scalar> &x )
{
	linear_combination(
		1,Teuchos::arrayArg<Scalar>(alpha)(),Teuchos::arrayArg<const Vector<Scalar>*>(&x)()
		,Teuchos::ScalarTraits<Scalar>::zero(),y
		);
}

template<class Scalar>
void TSFCore::Vp_StV( Vector<Scalar>* v_lhs, const Scalar& alpha, const Vector<Scalar>& v_rhs )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vp_StV(...), Error!");
#endif
  RTOpPack::TOpAXPY<Scalar> axpy_op(alpha);
	const Vector<Scalar>* vecs[]      = { &v_rhs };
	Vector<Scalar>*       targ_vecs[] = { v_lhs  };
	applyOp<Scalar>(axpy_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::Vp_V( Vector<Scalar>* y, const Vector<Scalar>& x, const Scalar& beta )
{
	linear_combination(
		1,Teuchos::arrayArg<Scalar>(Teuchos::ScalarTraits<Scalar>::one())()
		,Teuchos::arrayArg<const Vector<Scalar>*>(&x)()
		,beta,y
		);
}

template<class Scalar>
void TSFCore::abs( Vector<Scalar>* y, const Vector<Scalar>& x )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(y==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpAbs<Scalar> abs_op;
	const Vector<Scalar>* vecs[]      = { &x };
	Vector<Scalar>*       targ_vecs[] = { y  };
	applyOp<Scalar>(abs_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::reciprocal( Vector<Scalar>* y, const Vector<Scalar>& x )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(y==NULL,std::logic_error,"assign(...), Error!");
#endif
  RTOpPack::TOpReciprocal<Scalar> recip_op;
	const Vector<Scalar>* vecs[]      = { &x };
	Vector<Scalar>*       targ_vecs[] = { y  };
	applyOp<Scalar>(recip_op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
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
  RTOpPack::TOpEleWiseProd<Scalar> ele_wise_prod_op(alpha);
	const Vector<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
	Vector<Scalar>*       targ_vecs[] = { v_lhs };
	applyOp<Scalar>(ele_wise_prod_op,2,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
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
  RTOpPack::TOpEleWiseDivide<Scalar> ele_wise_divide_op(alpha);
	const Vector<Scalar>* vecs[]      = { &v_rhs1, &v_rhs2 };
	Vector<Scalar>*       targ_vecs[] = { v_lhs };
	applyOp<Scalar>(ele_wise_divide_op,2,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::linear_combination(
	const int                m
	,const Scalar            alpha[]
	,const Vector<Scalar>*   x[]
	,const Scalar            &beta
	,Vector<Scalar>          *y
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(y==NULL,std::logic_error,"linear_combination(...), Error!");
#endif
	if( beta == Teuchos::ScalarTraits<Scalar>::one() && m == 1 ) {
		Vp_StV( y, alpha[0], *x[0] );
		return;
	}
	else if( m == 0 ) {
		Vt_S( y, beta );
		return;
	}
  RTOpPack::TOpLinearCombination<Scalar> lin_comb_op(m,alpha,beta);
	Vector<Scalar>* targ_vecs[] = { y };
	applyOp<Scalar>(lin_comb_op,m,x,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

template<class Scalar>
void TSFCore::seed_randomize( unsigned int s )
{
  Teuchos::ScalarTraits<Scalar>::seedrandom(s);
}

template<class Scalar>
void TSFCore::randomize( Scalar l, Scalar u, Vector<Scalar>* v )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(v==NULL,std::logic_error,"Vt_S(...), Error");
#endif
  RTOpPack::TOpRandomize<Scalar> random_vector_op(l,u);
	Vector<Scalar>* targ_vecs[] = { v };
	applyOp<Scalar>(random_vector_op,0,(const Vector<Scalar>**)NULL,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
}

#endif // TSFCORE_VECTOR_STD_OPS_HPP



