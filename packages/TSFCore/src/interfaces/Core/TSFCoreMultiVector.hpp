// ///////////////////////////////////////////////////////////
// TSFCoreMultiVector.hpp

#ifndef TSFCORE_MULTI_VECTOR_HPP
#define TSFCORE_MULTI_VECTOR_HPP

#include "TSFCoreMultiVectorDecl.hpp"
//#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreAssertOp.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOpCppC.hpp"
#include "WorkspacePack.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

// Provide access to the columns as Vector objects

template<class Scalar>
Teuchos::RefCountPtr<const Vector<Scalar> >
MultiVector<Scalar>::col(Index j) const
{
	return const_cast<MultiVector*>(this)->col(j);
}

// Cloning

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
MultiVector<Scalar>::clone_mv() const
{
	return const_cast<MultiVector<Scalar>*>(this)->clone_mv();
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> > 
MultiVector<Scalar>::clone_mv()
{
	const VectorSpace<Scalar>
		&domain = *this->domain(),
		&range  = *this->range();
	Teuchos::RefCountPtr<MultiVector<Scalar> >
		copy = range.createMembers(domain.dim());
	RTOpPack::RTOpC assign_vectors_op;
	if(0>RTOp_TOp_assign_vectors_construct(&assign_vectors_op.op())) assert(0);
	const MultiVector<Scalar>* multi_vecs[1]      = { this };
	MultiVector<Scalar>*       targ_multi_vecs[1] = { copy.get() };
	TSFCore::applyOp<Scalar>(assign_vectors_op,1,multi_vecs,1,targ_multi_vecs,RTOp_REDUCT_OBJ_NULL);
	return copy;
}

// Sub-view methods

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
MultiVector<Scalar>::subView( const Range1D& col_rng ) const
{
	return const_cast<MultiVector<Scalar>*>(this)->subView(col_rng);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MultiVector<Scalar>::subView( const Range1D& col_rng_in )
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore        *wss      = WorkspacePack::default_workspace_store.get();
	const VectorSpace<Scalar>  &domain   = *this->domain();
	const VectorSpace<Scalar>  &range    = *this->range();
	const Index                dimDomain = domain.dim();
	const Range1D              col_rng   = RangePack::full_range(col_rng_in,1,dimDomain);
	if( col_rng.lbound() == 1 && col_rng.ubound() == dimDomain )
		return Teuchos::rcp(this,false); // Takes all of the colunns!
	if( col_rng.size() ) {
		// We have to create a view of a subset of the columns
		wsp::Workspace< Teuchos::RefCountPtr< Vector<Scalar> > >  col_vecs(wss,col_rng.size());
		for( Index j = col_rng.lbound(); j <= col_rng.ubound(); ++j )
			col_vecs[j-1] = this->col(j);
		return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(col_rng.size()),&col_vecs[0]));
	}
	return Teuchos::null; // There was an empty set in col_rng_in!
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
MultiVector<Scalar>::subView( const int numCols, const int cols[] ) const
{
	return const_cast<MultiVector<Scalar>*>(this)->subView(numCols,cols);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MultiVector<Scalar>::subView( const int numCols, const int cols[] )
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore        *wss      = WorkspacePack::default_workspace_store.get();
	const VectorSpace<Scalar>  &domain   = *this->domain();
	const VectorSpace<Scalar>  &range    = *this->range();
	const Index                dimDomain = domain.dim();
#ifdef _DEBUG
	const char msg_err[] = "MultiVector<Scalar>::subView(numCols,cols[]): Error!";
 	TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
	// We have to create a view of a subset of the columns
	wsp::Workspace< Teuchos::RefCountPtr< Vector<Scalar> > > col_vecs(wss,numCols);
	for( int k = 0; k < numCols; ++k ) {
		const int col_k = cols[k];
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			col_k < 1 || dimDomain < col_k, std::invalid_argument
			,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [1,"<<dimDomain<<"]!"
			);
#endif
		col_vecs[k] = this->col(col_k);
	}
	return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(numCols),&col_vecs[0]));
}

// Collective applyOp() methods

template<class Scalar>
void MultiVector<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &prim_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_objs[]
	,const Index                    prim_first_ele_in
	,const Index                    prim_sub_dim_in
	,const Index                    prim_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	const VectorSpace<Scalar>
		&domain = *this->domain(),
		&range  = *this->range();

	// Get the primary and secondary dimmensions.
	const Index
		prim_dim     = range.dim(),
		sec_dim      = domain.dim(),
		prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim - prim_first_ele_in + 1 ),
		sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  -  sec_first_ele_in + 1  );
#ifdef _DEBUG
	const char err_msg[] = "MultiVector<Scalar>::applyOp(...): Error!";
	TEST_FOR_EXCEPTION( !(0 < prim_sub_dim && prim_sub_dim <= prim_dim), std::invalid_argument, err_msg );
	TEST_FOR_EXCEPTION( !(0 < sec_sub_dim  && sec_sub_dim  <= sec_dim),  std::invalid_argument, err_msg );
#endif

	//
	// Apply the reduction/transformation operator and transform the
	// target vectors and reduce each of the reduction objects.
	//

	wsp::Workspace< Teuchos::RefCountPtr<const Vector<Scalar> > >   vecs_s(wss,num_multi_vecs);
	wsp::Workspace<const Vector<Scalar>*>                         vecs(wss,num_multi_vecs);
	wsp::Workspace< Teuchos::RefCountPtr<Vector<Scalar> > >         targ_vecs_s(wss,num_targ_multi_vecs);
	wsp::Workspace<Vector<Scalar>*>                               targ_vecs(wss,num_targ_multi_vecs);

	for(Index j = sec_first_ele_in; j <= sec_first_ele_in - 1 + sec_sub_dim; ++j) {
		// Fill the arrays of vector arguments 
		{for(Index k = 0; k < num_multi_vecs; ++k) {
			vecs_s[k] = multi_vecs[k]->col(j);
			vecs[k] = vecs_s[k].get();
		}}
		{for(Index k = 0; k < num_targ_multi_vecs; ++k) {
			targ_vecs_s[k] = targ_multi_vecs[k]->col(j);
			targ_vecs[k] = targ_vecs_s[k].get();
		}}
		// Apply the reduction/transformation operator
		TSFCore::applyOp(
			prim_op
			,num_multi_vecs,      (num_multi_vecs      ? &vecs[0]      : NULL)
			,num_targ_multi_vecs, (num_targ_multi_vecs ? &targ_vecs[0] : NULL)
			,reduct_objs ? reduct_objs[j-1] : NULL
			,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
			);
	}
	// At this point all of the designated targ vectors in the target multi-vectors have
	// been transformed and all the reduction objects in reduct_obj[] have accumulated
	// the reductions.
}

template<class Scalar>
void MultiVector<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &prim_op
	,const RTOpPack::RTOpT<Scalar>  &sec_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_obj
	,const Index                    prim_first_ele_in
	,const Index                    prim_sub_dim_in
	,const Index                    prim_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	const VectorSpace<Scalar>
		&domain = *this->domain(),
		&range  = *this->range();

	// Get the primary and secondary dimmensions.
	const Index
		prim_dim     = range.dim(),
		sec_dim      = domain.dim(),
		prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim - prim_first_ele_in + 1 ),
		sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  -  sec_first_ele_in + 1  );
#ifdef _DEBUG
	const char err_msg[] = "MultiVector<Scalar>::applyOp(...): Error!";
	TEST_FOR_EXCEPTION( !(0 < prim_sub_dim && prim_sub_dim <= prim_dim), std::invalid_argument, err_msg );
	TEST_FOR_EXCEPTION( !(0 < sec_sub_dim  && sec_sub_dim  <= sec_dim),  std::invalid_argument, err_msg );
#endif

	// Create a temporary buffer for the reduction objects of the primary reduction
	// so that we can call the companion version of this method.
	wsp::Workspace<RTOp_ReductTarget>   reduct_objs(wss,reduct_obj!=RTOp_REDUCT_OBJ_NULL?sec_sub_dim:0);
	if(reduct_obj!=RTOp_REDUCT_OBJ_NULL) {
		for(Index k = 0; k < sec_sub_dim; ++k) {
			prim_op.reduct_obj_create_raw( &(reduct_objs[k]=RTOp_REDUCT_OBJ_NULL) );
		}
	}
	
	// Call the campanion version that accepts an array of reduction objects
	this->applyOp(
		prim_op
		,num_multi_vecs,       multi_vecs
		,num_targ_multi_vecs,  targ_multi_vecs
		,reduct_obj!=RTOp_REDUCT_OBJ_NULL ? &reduct_objs[0] : NULL
		,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
		,sec_first_ele_in,  sec_sub_dim_in
		);

	// Reduce all the reduction objects using the secondary reduction operator
	// into one reduction object and free the intermedate reduction objects.
	if(reduct_obj!=RTOp_REDUCT_OBJ_NULL) {
		for(Index k = 0; k < sec_sub_dim; ++k) {
			sec_op.reduce_reduct_objs( reduct_objs[k] ,reduct_obj );
			prim_op.reduct_obj_free( &reduct_objs[k] );
		}
	}
}

// Overridden methods from LinearOp

template<class Scalar>
void MultiVector<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	Vp_MtV_assert_compatibility(y,*this,M_trans,x);
	const Index nc = this->domain()->dim();
	// y *= beta
	Vt_S(y,beta);
	// y += alpha*op(M)*x
	if(M_trans == NOTRANS) {
		//
		// y += alpha*M*x = alpha*M.col(1)*x(1) + ... + alpha*M.col(nc)*x(nc)
		//
		// Extract an explicit view of x
		RTOpPack::SubVectorT<Scalar> x_sub_vec;               
		x.getSubVector(Range1D(),&x_sub_vec);
		// Loop through and add the multiple if each column
		for(Index j = 1; j <= nc; ++j )
			Vp_StV( y, alpha*x_sub_vec(j), *this->col(j) );
		// Release the view of x
		x.freeSubVector(&x_sub_vec);
	}
	else {
		//
		//                   [ alpha*dot(M.col(1),x)  ]
		// y += alpha*M'*x = [ alpha*dot(M.col(2),x)  ]
		//                   [ ...                    ]
		//                   [ alpha*dot(M.col(nc),x) ]
		//
		// Extract an explicit view of y
		RTOpPack::MutableSubVectorT<Scalar> y_sub_vec;               
		y->getSubVector(Range1D(),&y_sub_vec);
		// Loop through and add to each element in y
		for(Index j = 1; j <= nc; ++j )
			y_sub_vec(j) += alpha*dot(*this->col(j),x);
		// Commit explicit view of y
		y->commitSubVector(&y_sub_vec);
	}
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
MultiVector<Scalar>::clone() const
{
	return this->clone_mv();
}

} // end namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_HPP
