// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIMultiVectorBase.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_HPP

#include <vector>

#include "TSFCoreMPIMultiVectorBaseDecl.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"
#include "RTOpCppToMPI.hpp"
#include "WorkspacePack.hpp"
#include "dynamic_cast_verbose.hpp"

namespace TSFCore {

template<class Scalar>
MPIMultiVectorBase<Scalar>::MPIMultiVectorBase()
	:in_applyOp_(false)
	,globalDim_(0)
	,localOffset_(-1)
	,localSubDim_(0)
	,numCols_(0)
{}

// Overridden form OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
MPIMultiVectorBase<Scalar>::range() const
{
	return mpiSpace();
}

// Overridden from LinearOp

template<class Scalar>
void MPIMultiVectorBase<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	MultiVectorCols<Scalar>
		X(Teuchos::rcp(const_cast<Vector<Scalar>*>(&x),false)),
		Y(Teuchos::rcp(y,false));
	apply(M_trans,X,&Y,alpha,beta);
}

// Overridden from MultiVector

template<class Scalar>
void MPIMultiVectorBase<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &pri_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_objs[]
	,const Index                    pri_first_ele_in
	,const Index                    pri_sub_dim_in
	,const Index                    pri_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	const Index numCols = this->domain()->dim();
	const MPIVectorSpaceBase<Scalar> &mpiSpc = *mpiSpace();
#ifdef _DEBUG
	// ToDo: Validate input!
	TEST_FOR_EXCEPTION(
		in_applyOp_, std::invalid_argument
		,"MPIMultiVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
		"clear sign that one of the methods getSubMultiVector(...), freeSubMultiVector(...) or commitSubMultiVector(...) "
		"was not implemented properly!"
		);
#endif
	// Flag that we are in applyOp()
	in_applyOp_ = true;
	// Get the overlap in the current process with the input logical sub-vector
	// from (first_ele_in,sub_dim_in,global_offset_in)
	RTOp_index_type  overlap_first_local_ele  = 0;
	RTOp_index_type  overalap_local_sub_dim   = 0;
	RTOp_index_type  overlap_global_offset    = 0;
	RTOp_parallel_calc_overlap(
		globalDim_, localSubDim_, localOffset_, pri_first_ele_in, pri_sub_dim_in, pri_global_offset_in
		,&overlap_first_local_ele, &overalap_local_sub_dim, &overlap_global_offset
		);
	const Range1D
		local_rng = (
			overlap_first_local_ele!=0
			? Range1D( localOffset_ + overlap_first_local_ele, localOffset_ + overlap_first_local_ele + overalap_local_sub_dim - 1 )
			: Range1D::Invalid
			),
		col_rng(
			sec_first_ele_in
			,sec_sub_dim_in ? sec_first_ele_in + sec_sub_dim_in - 1 : numCols
			);
	// Create sub-vector views of all of the *participating* local data
	wsp::Workspace<RTOpPack::SubMultiVectorT<Scalar> > sub_multi_vecs(wss,num_multi_vecs);
	wsp::Workspace<RTOpPack::MutableSubMultiVectorT<Scalar> > targ_sub_multi_vecs(wss,num_targ_multi_vecs);
	if( overlap_first_local_ele != 0 ) {
		if(1){for(int k = 0; k < num_multi_vecs; ++k ) {
			multi_vecs[k]->getSubMultiVector( local_rng, col_rng, &sub_multi_vecs[k] );
			sub_multi_vecs[k].setGlobalOffset( overlap_global_offset );
		}}
		if(1){for(int k = 0; k < num_targ_multi_vecs; ++k ) {
			targ_multi_vecs[k]->getSubMultiVector( local_rng, col_rng, &targ_sub_multi_vecs[k] );
			targ_sub_multi_vecs[k].setGlobalOffset( overlap_global_offset );
		}}
	}
	// Apply the RTOp operator object (all processors must participate)
	RTOpPack::MPI_apply_op(
		mpiSpc.mpiComm()                                                                   // comm
		,pri_op                                                                            // op
		,-1                                                                                // root_rank (perform an all-reduce)
		,col_rng.size()                                                                    // num_cols
		,num_multi_vecs                                                                    // num_multi_vecs
		,num_multi_vecs && overlap_first_local_ele ? &sub_multi_vecs[0] : NULL             // sub_multi_vecs
		,num_targ_multi_vecs                                                               // num_targ_multi_vecs
		,num_targ_multi_vecs && overlap_first_local_ele ? &targ_sub_multi_vecs[0] : NULL   // targ_sub_multi_vecs
		,reduct_objs                                                                       // reduct_objs
		);
	// Free and commit the local data
	if(1){for(int k = 0; k < num_multi_vecs; ++k ) {
		sub_multi_vecs[k].setGlobalOffset(local_rng.lbound()-1);
		multi_vecs[k]->freeSubMultiVector( &sub_multi_vecs[k] );
	}}
	if(1){for(int k = 0; k < num_targ_multi_vecs; ++k ) {
		targ_sub_multi_vecs[k].setGlobalOffset(local_rng.lbound()-1);
		targ_multi_vecs[k]->commitSubMultiVector( &targ_sub_multi_vecs[k] );
	}}
	// Flag that we are leaving applyOp()
	in_applyOp_ = false;
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::getSubMultiVector(
	const Range1D                       &rowRng_in
	,const Range1D                      &colRng_in
	,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
	) const
{
	const Range1D rowRng = validateRowRange(rowRng_in);
	const Range1D colRng = validateColRange(colRng_in);
	if( rowRng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rowRng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		MultiVector<Scalar>::getSubMultiVector(rowRng_in,colRng_in,sub_mv);
		return;
	}
	// rng consists of all local data so get it!
	const Scalar *localValues = NULL;
	int leadingDim = 0;
	this->getLocalData(&localValues,&leadingDim);
	sub_mv->initialize(
		rowRng.lbound()-1                             // globalOffset
		,rowRng.size()                                // subDim
		,colRng.lbound()-1                            // colOffset
		,colRng.size()                                // numSubCols
		,localValues
		+(rowRng.lbound()-localOffset_-1)
		+(colRng.lbound()-1)*leadingDim               // values
		,leadingDim                                   // leadingDim
		);
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::freeSubMultiVector(
	RTOpPack::SubMultiVectorT<Scalar>* sub_mv
	) const
{
	if( sub_mv->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_mv->globalOffset()+sub_mv->subDim() ) {
		// Let the default implementation handle it!
		MultiVector<Scalar>::freeSubMultiVector(sub_mv);
		return;
	}
	freeLocalData( sub_mv->values() );
	sub_mv->set_uninitialized();
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::getSubMultiVector(
	const Range1D                                &rowRng_in
	,const Range1D                               &colRng_in
	,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
	)
{
	const Range1D rowRng = validateRowRange(rowRng_in);
	const Range1D colRng = validateColRange(colRng_in);
	if( rowRng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rowRng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		MultiVector<Scalar>::getSubMultiVector(rowRng_in,colRng_in,sub_mv);
		return;
	}
	// rng consists of all local data so get it!
	Scalar *localValues = NULL;
	int leadingDim = 0;
	this->getLocalData(&localValues,&leadingDim);
	sub_mv->initialize(
		rowRng.lbound()-1                             // globalOffset
		,rowRng.size()                                // subDim
		,colRng.lbound()-1                            // colOffset
		,colRng.size()                                // numSubCols
		,localValues
		+(rowRng.lbound()-localOffset_-1)
		+(colRng.lbound()-1)*leadingDim               // values
		,leadingDim                                   // leadingDim
		);
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::commitSubMultiVector(
	RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv
	)
{
	if( sub_mv->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_mv->globalOffset()+sub_mv->subDim() ) {
		// Let the default implementation handle it!
		MultiVector<Scalar>::commitSubMultiVector(sub_mv);
		return;
	}
	commitLocalData( sub_mv->values() );
	sub_mv->set_uninitialized();
}

// protected

template<class Scalar>
void MPIMultiVectorBase<Scalar>::updateMpiSpace()
{
	if(globalDim_ == 0) {
		const MPIVectorSpaceBase<Scalar> *mpiSpace = this->mpiSpace().get();
		if(mpiSpace) {
			globalDim_    = mpiSpace->dim();
			localOffset_  = mpiSpace->localOffset();
			localSubDim_  = mpiSpace->localSubDim();
			numCols_      = this->domain()->dim();
		}
		else {
			globalDim_    = 0;
			localOffset_  = -1;
			localSubDim_  = 0;
			numCols_      = 0;
		}
	}
}

template<class Scalar>
Range1D MPIMultiVectorBase<Scalar>::validateRowRange( const Range1D &rowRng_in ) const
{
	const Range1D rowRng = RangePack::full_range(rowRng_in,1,globalDim_);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		rowRng.lbound() < 1 || globalDim_ < rowRng.ubound(), std::invalid_argument
		,"MPIMultiVectorBase<Scalar>::validateRowRange(rowRng): Error, the range rowRng = ["
		<<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not "
		"in the range [1,"<<globalDim_<<"]!"
		);
#endif
	return rowRng;
}

template<class Scalar>
Range1D MPIMultiVectorBase<Scalar>::validateColRange( const Range1D &colRng_in ) const
{
	const Range1D colRng = RangePack::full_range(colRng_in,1,numCols_);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		colRng.lbound() < 1 || numCols_ < colRng.ubound(), std::invalid_argument
		,"MPIMultiVectorBase<Scalar>::validateColRange(colRng): Error, the range colRng = ["
		<<colRng.lbound()<<","<<colRng.ubound()<<"] is not "
		"in the range [1,"<<numCols_<<"]!"
		);
#endif
	return colRng;
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_HPP
