// ////////////////////////////////////////////////////////////////////
// TSFCoreMultiVectorStdOps.hpp

#ifndef TSFCORE_MULTI_VECTOR_STD_OPS_HPP
#define TSFCORE_MULTI_VECTOR_STD_OPS_HPP

#include "TSFCoreMultiVectorStdOpsDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_ROp_max.h"
#include "RTOp_ROp_sum_abs.h"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOpCppC.hpp"
#include "Teuchos_TestForException.hpp"

template<class Scalar>
void TSFCore::dot( const MultiVector<Scalar>& V1, const MultiVector<Scalar>& V2, Scalar dot[] )
{
	int kc;
    const int m = V1.domain()->dim();
    RTOpPack::RTOpC dot_op;
    RTOp_ROp_dot_prod_construct(&dot_op.op());
    std::vector<RTOp_ReductTarget>  dot_targs(m);
    for( kc = 0; kc < m; ++kc )
        dot_op.reduct_obj_create_raw(&(dot_targs[kc]=RTOp_REDUCT_OBJ_NULL));
    const MultiVector<Scalar>* multi_vecs[] = { &V1, &V2 };
    applyOp<Scalar>(dot_op,2,multi_vecs,0,NULL,&dot_targs[0]);
    for( kc = 0; kc < m; ++kc ) {
        dot[kc] = RTOp_ROp_dot_prod_val(dot_targs[kc]);
        dot_op.reduct_obj_free(&(dot_targs[kc]));
    }
}

template<class Scalar>
Scalar TSFCore::norm_1( const MultiVector<Scalar>& V )
{
	// Primary column-wise reduction (sum of absolute values)
	RTOpPack::RTOpC sum_abs_op;
	RTOp_ROp_sum_construct(&sum_abs_op.op());
	// Secondaary reduction (max over all columns = norm_1)
	RTOpPack::RTOpC max_op;
	RTOp_ROp_max_construct(&max_op.op());
	// Reduction object (must be same for both sum_abs and max_targ objects)
	RTOpPack::ReductTarget max_targ;
	max_op.reduct_obj_create(&max_targ);
	// Perform the reductions
    const MultiVector<Scalar>* multi_vecs[] = { &V };
    applyOp<Scalar>(sum_abs_op,max_op,1,multi_vecs,0,NULL,max_targ.obj());
	// Return the final value
	return RTOp_ROp_max_val(max_targ.obj());
}

template<class Scalar>
void TSFCore::scale( Scalar alpha, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::RTOpC  scale_vector_op;
	if(0>RTOp_TOp_scale_vector_construct(alpha,&scale_vector_op.op())) assert(0);
	MultiVector<Scalar>* targ_multi_vecs[] = { V };
	applyOp<Scalar>(scale_vector_op,0,NULL,1,targ_multi_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::assign( MultiVector<Scalar>* V, Scalar alpha )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::RTOpC  assign_scalar_op;
	if(0>RTOp_TOp_assign_scalar_construct(alpha,&assign_scalar_op.op())) assert(0);
	MultiVector<Scalar>* targ_multi_vecs[] = { V };
	applyOp<Scalar>(assign_scalar_op,0,NULL,1,targ_multi_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::assign( MultiVector<Scalar>* V, const MultiVector<Scalar>& U )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"assign(...), Error!");
#endif
	RTOpPack::RTOpC assign_vectors_op;
	if(0>RTOp_TOp_assign_vectors_construct(&assign_vectors_op.op())) assert(0);
	const MultiVector<Scalar>* multi_vecs[]      = { &U };
	MultiVector<Scalar>*       targ_multi_vecs[] = { V   };
	applyOp<Scalar>(assign_vectors_op,1,multi_vecs,1,targ_multi_vecs,RTOp_REDUCT_OBJ_NULL);
}

template<class Scalar>
void TSFCore::update( Scalar alpha, const MultiVector<Scalar>& U, MultiVector<Scalar>* V )
{
#ifdef _DEBUG
    TEST_FOR_EXCEPTION(V==NULL,std::logic_error,"update(...), Error!");
#endif
    RTOpPack::RTOpC  axpy_op;
    RTOp_TOp_axpy_construct(alpha,&axpy_op.op());
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
