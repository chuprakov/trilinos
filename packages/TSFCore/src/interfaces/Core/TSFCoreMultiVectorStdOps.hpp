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
// TSFCoreMultiVectorStdOps.hpp

#ifndef TSFCORE_MULTI_VECTOR_STD_OPS_HPP
#define TSFCORE_MULTI_VECTOR_STD_OPS_HPP

#include "TSFCoreMultiVectorStdOpsDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "RTOpPack_TOpAssignVectors.hpp"
#include "RTOpPack_TOpAXPY.hpp"
#include "RTOpPack_TOpScaleVector.hpp"
#include "Teuchos_TestForException.hpp"

template<class Scalar>
void TSFCore::norms( const MultiVector<Scalar>& V, Scalar norms[] )
{
	const int m = V.domain()->dim();
	V.range()->scalarProds(V,V,norms);
	for( int j = 0; j < m; ++j )
		norms[j] = Teuchos::ScalarTraits<Scalar>::squareroot(norms[j]);
}

template<class Scalar, class NormOp>
void TSFCore::reductions( const MultiVector<Scalar>& V, const NormOp &op, Scalar norms[] )
{
	int kc;
  const int m = V.domain()->dim();
  std::vector<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >  rcp_op_targs(m);
  std::vector<RTOpPack::ReductTarget*>                        op_targs(m);
  for( kc = 0; kc < m; ++kc ) {
    rcp_op_targs[kc] = op.reduct_obj_create();
    op_targs[kc] = &*rcp_op_targs[kc];
  }
  const MultiVector<Scalar>* multi_vecs[] = { &V,};
  applyOp<Scalar>(op,1,multi_vecs,0,NULL,&op_targs[0]);
  for( kc = 0; kc < m; ++kc ) {
    norms[kc] = op(*op_targs[kc]);
  }
}

template<class Scalar>
void TSFCore::dots( const MultiVector<Scalar>& V1, const MultiVector<Scalar>& V2, Scalar dots[] )
{
	int kc;
  const int m = V1.domain()->dim();
  RTOpPack::ROpDotProd<Scalar> dot_op;
  std::vector<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >  rcp_dot_targs(m);
  std::vector<RTOpPack::ReductTarget*>                        dot_targs(m);
  for( kc = 0; kc < m; ++kc ) {
    rcp_dot_targs[kc] = dot_op.reduct_obj_create();
    dot_targs[kc] = &*rcp_dot_targs[kc];
  }
  const MultiVector<Scalar>* multi_vecs[] = { &V1, &V2 };
  applyOp<Scalar>(dot_op,2,multi_vecs,0,NULL,&dot_targs[0]);
  for( kc = 0; kc < m; ++kc ) {
    dots[kc] = dot_op(*dot_targs[kc]);
  }
}

template<class Scalar>
Scalar TSFCore::norm_1( const MultiVector<Scalar>& V )
{
	// Primary column-wise reduction (sum of absolute values)
  RTOpPack::ROpNorm1<Scalar> sum_abs_op;
	// Secondaary reduction (max over all columns = induced norm_1 matrix  norm)
  RTOpPack::ROpNormInf<Scalar> max_op;
	// Reduction object (must be same for both sum_abs and max_targ objects)
  Teuchos::RefCountPtr<RTOpPack::ReductTarget>
    max_targ = max_op.reduct_obj_create();
	// Perform the reductions
  const MultiVector<Scalar>* multi_vecs[] = { &V };
  applyOp<Scalar>(sum_abs_op,max_op,1,multi_vecs,0,NULL,&*max_targ);
	// Return the final value
	return max_op(*max_targ);
}

template<class Scalar>
void TSFCore::scale( Scalar alpha, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
	if(alpha==Teuchos::ScalarTraits<Scalar>::zero()) {
		assign( V, Teuchos::ScalarTraits<Scalar>::zero() );
		return;
	}
	if(alpha==Teuchos::ScalarTraits<Scalar>::one()) {
		return;
	}
  RTOpPack::TOpScaleVector<Scalar> scale_vector_op(alpha);
	MultiVector<Scalar>* targ_multi_vecs[] = { V };
	applyOp<Scalar>(
		scale_vector_op,0,(const MultiVector<Scalar>**)NULL // The SUN compiler requires these casts!
		,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL
		);
}

template<class Scalar>
void TSFCore::scaleUpdate( const Vector<Scalar>& a, const MultiVector<Scalar>& U, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
	bool is_compatible = U.range()->isCompatible(*a.space());
	TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.range()->isCompatible(*a.space())==false");
	is_compatible = U.range()->isCompatible(*V->range());
	TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.range()->isCompatible((V->range())==false ");
	is_compatible = U.domain()->isCompatible(*V->domain());
	TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.domain().isCompatible(V->domain())==false ");
#endif
	const int m = U.domain()->dim();
	for( int j = 1; j <= m; ++j ) {
		ele_wise_prod( Scalar(1.0), a, *U.col(j), &*V->col(j) ); 
	}
}

template<class Scalar>
void TSFCore::assign( MultiVector<Scalar>* V, Scalar alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
	MultiVector<Scalar>* targ_multi_vecs[] = { V };
	applyOp<Scalar>(
		assign_scalar_op,0,(const MultiVector<Scalar>**)NULL // The SUN compiler requires these casts!
		,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL
		);
}

template<class Scalar>
void TSFCore::assign( MultiVector<Scalar>* V, const MultiVector<Scalar>& U )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::TOpAssignVectors<Scalar> assign_vectors_op;
	const MultiVector<Scalar>* multi_vecs[]      = { &U };
	MultiVector<Scalar>*       targ_multi_vecs[] = { V   };
	applyOp<Scalar>(
		assign_vectors_op,1,multi_vecs,1,targ_multi_vecs
		,(RTOpPack::ReductTarget**)NULL // The SUN compiler requires this cast!
		);
}

template<class Scalar>
void TSFCore::update( Scalar alpha, const MultiVector<Scalar>& U, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
#endif
  RTOpPack::TOpAXPY<Scalar> axpy_op(alpha);
  const MultiVector<Scalar>* multi_vecs[]       = { &U };
  MultiVector<Scalar>*       targ_multi_vecs[]  = { V  };
  applyOp<Scalar>(axpy_op,1,multi_vecs,1,targ_multi_vecs,NULL);
}

template<class Scalar>
void TSFCore::update( Scalar alpha[], Scalar beta, const MultiVector<Scalar>& U, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
	bool is_compatible = U.range()->isCompatible(*V->range());
  TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.range()->isCompatible((V->range())==false ");
	is_compatible = U.domain()->isCompatible(*V->domain());
  TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.domain().isCompatible(V->domain())==false ");
#endif
	const int m = U.domain()->dim();
	for( int j = 1; j <= m; ++j )
		Vp_StV( V->col(j).get(), alpha[j-1]*beta, *U.col(j) );
}


template<class Scalar>
void TSFCore::update( const MultiVector<Scalar>& U, Scalar alpha[], Scalar beta, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
    TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
	bool is_compatible = U.range()->isCompatible(*V->range());
    TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.range()->isCompatible((V->range())==false ");
	is_compatible = U.domain()->isCompatible(*V->domain());
    TEST_FOR_EXCEPTION(
		!is_compatible,Exceptions::IncompatibleVectorSpaces
		,"update(...), Error, U.domain().isCompatible(V->domain())==false ");
#endif
	const int m = U.domain()->dim();
	for( int j = 1; j <= m; ++j ) {
		Vt_S( V->col(j).get(), alpha[j-1]*beta );
		Vp_StV( V->col(j).get(), 1.0, *U.col(j) );
	}
}

template<class Scalar>
void TSFCore::linear_combination(
	const int                     m
	,const Scalar                 alpha[]
	,const MultiVector<Scalar>*   X[]
	,const Scalar                 &beta
	,MultiVector<Scalar>          *Y
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(Y==NULL,std::logic_error,"linear_combination(...), Error!");
#endif
	if( beta == Teuchos::ScalarTraits<Scalar>::one() && m == 1 ) {
		update( alpha[0], *X[0], Y );
		return;
	}
	else if( m == 0 ) {
		scale( beta, Y );
		return;
	}
  RTOpPack::TOpLinearCombination<Scalar> lin_comb_op(m,alpha,beta);
	MultiVector<Scalar>* targ_multi_vecs[] = { Y };
	TSFCore::applyOp<Scalar>(lin_comb_op,m,X,1,targ_multi_vecs,(RTOpPack::ReductTarget**)NULL);
}

template<class Scalar>
void TSFCore::randomize( Scalar l, Scalar u, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
    TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"randomize(...), Error!");
#endif
	const int m = V->domain()->dim();
	for( int j = 1; j <= m; ++j ) {
		randomize( l, u, V->col(j).get() ); // Todo: call applyOp(...) directly!
	}
}

#endif // TSFCORE_MULTI_VECTOR_STD_OPS_HPP
