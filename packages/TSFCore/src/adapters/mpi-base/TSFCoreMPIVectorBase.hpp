// /////////////////////////////////////////////////////////////////
// TSFCoreMPIVectorBase.hpp

#ifndef TSFCORE_MPI_VECTOR_BASE_HPP
#define TSFCORE_MPI_VECTOR_BASE_HPP

#include <stdexcept>

#include "TSFCoreMPIVectorBaseDecl.hpp"
#include "RTOp_parallel_helpers.h"
#include "RTOpCppToMPI.hpp"
#include "WorkspacePack.hpp"
#include "ThrowException.hpp"
#include "dynamic_cast_verbose.hpp"

namespace TSFCore {

template<class Scalar>
MPIVectorBase<Scalar>::MPIVectorBase()
	:in_applyOp_(false)
	,globalDim_(-1)
	,localOffset_(-1)
	,localSubDim_(-1)
{}

// Virtual methods with default implementations

template<class Scalar>
void MPIVectorBase<Scalar>::getLocalData( const Scalar** values, ptrdiff_t* stride ) const
{
	const_cast<MPIVectorBase<Scalar>*>(this)->getLocalData((Scalar**)values,stride);
}

// Overridden from Vector

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >
MPIVectorBase<Scalar>::space() const
{
	return mpiSpace();
}

template<class Scalar>
void MPIVectorBase<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &op
	,const size_t                   num_vecs
	,const Vector<Scalar>*          vecs[]
	,const size_t                   num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOp_ReductTarget              reduct_obj
	,const Index                    first_ele_in
	,const Index                    sub_dim_in
	,const Index                    global_offset_in
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	namespace mmp = MemMngPack;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	const MPIVectorSpaceBase<Scalar>  &mpiSpc = *mpiSpace();
	const Index globalDim = mpiSpc.dim(), localOffset = mpiSpc.localOffset(), localSubDim = mpiSpc.localSubDim();
#ifdef _DEBUG
	// ToDo: Validate input!
	THROW_EXCEPTION(
		in_applyOp_, std::invalid_argument
		,"MPIVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
		"clear sign that one of the methods getSubVector(...), freeSubVector(...) or commitSubVector(...) "
		"was not implemented properly!"
		);
#endif
	// Flag that we are in applyOp()
	in_applyOp_ = true;
	// Convert Vector<> arguments into MPIVectorBase<> objects
	wsp::Workspace<const MPIVectorBase<Scalar>*> mpi_vecs(wss,num_vecs);
	if(1){for(int k = 0; k < num_vecs; ++k )
		mpi_vecs[k] = &dyn_cast<const MPIVectorBase<Scalar> >(*vecs[k]);
	}
	wsp::Workspace<MPIVectorBase<Scalar>*> mpi_targ_vecs(wss,num_targ_vecs);
	if(1){for(int k = 0; k < num_targ_vecs; ++k )
		mpi_targ_vecs[k] = &dyn_cast<MPIVectorBase<Scalar> >(*targ_vecs[k]);
	}
	// Get the overlap in the current process with the input logical sub-vector
	// from (first_ele_in,sub_dim_in,global_offset_in)
	RTOp_index_type  overlap_first_local_ele  = 0;
	RTOp_index_type  overalap_local_sub_dim   = 0;
	RTOp_index_type  overlap_global_offset    = 0;
	RTOp_parallel_calc_overlap(
		globalDim, localSubDim, localOffset, first_ele_in, sub_dim_in, global_offset_in
		,&overlap_first_local_ele, &overalap_local_sub_dim, &overlap_global_offset
		);
	const Range1D local_rng = (
		overlap_first_local_ele!=0
		? Range1D( localOffset + overlap_first_local_ele, localOffset + overlap_first_local_ele + overalap_local_sub_dim - 1 )
		: Range1D::Invalid
		);
	// Create sub-vector views of all of the *participating* local data
	wsp::Workspace<RTOpPack::SubVectorT<Scalar> > sub_vecs(wss,num_vecs);
	wsp::Workspace<RTOpPack::MutableSubVectorT<Scalar> > sub_targ_vecs(wss,num_targ_vecs);
	if( overlap_first_local_ele != 0 ) {
		if(1){for(int k = 0; k < num_vecs; ++k ) {
			mpi_vecs[k]->getSubVector( local_rng, &sub_vecs[k] );
			sub_vecs[k].setGlobalOffset( overlap_global_offset );
		}}
		if(1){for(int k = 0; k < num_targ_vecs; ++k ) {
			mpi_targ_vecs[k]->getSubVector( local_rng, &sub_targ_vecs[k] );
			sub_targ_vecs[k].setGlobalOffset( overlap_global_offset );
		}}
	}
	// Apply the RTOp operator object
	RTOp_ReductTarget reduct_objs[] = { reduct_obj };
	RTOpPack::MPI_apply_op(
		mpiSpc.mpiComm()                                                       // comm
		,op                                                                    // op
		,-1                                                                    // root_rank (perform an all-reduce)
		,1                                                                     // num_cols
		,num_vecs                                                              // num_vecs
		,num_vecs && overlap_first_local_ele ? &sub_vecs[0] : NULL             // sub_vecs
		,num_targ_vecs                                                         // num_targ_vecs
		,num_targ_vecs && overlap_first_local_ele ? &sub_targ_vecs[0] : NULL   // targ_sub_vecs
		,reduct_objs                                                           // reduct_objs
		);
	// Free and commit the local data
	if(1){for(int k = 0; k < num_vecs; ++k ) {
		sub_vecs[k].setGlobalOffset(local_rng.lbound()-1);
		mpi_vecs[k]->freeSubVector( &sub_vecs[k] );
	}}
	if(1){for(int k = 0; k < num_targ_vecs; ++k ) {
		sub_targ_vecs[k].setGlobalOffset(local_rng.lbound()-1);
		mpi_targ_vecs[k]->commitSubVector( &sub_targ_vecs[k] );
	}}
	// Flag that we are leaving applyOp()
	in_applyOp_ = false;
}

template<class Scalar>
void MPIVectorBase<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	const Range1D rng = validateRange(rng_in);
	if( rng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		Vector<Scalar>::getSubVector(rng_in,sub_vec);
		return;
	}
	// rng consists of all local data so get it!
	const Scalar *localValues = NULL;
	ptrdiff_t stride = 0;
	this->getLocalData(&localValues,&stride);
	sub_vec->initialize(
		rng.lbound()-1                             // globalOffset
		,rng.size()                                // subDim
		,localValues+(rng.lbound()-localOffset_-1) // values
		,stride                                    // stride
		);
}

template<class Scalar>
void MPIVectorBase<Scalar>::freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	if( sub_vec->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim() ) {
		// Let the default implementation handle it!
		Vector<Scalar>::freeSubVector(sub_vec);
		return;
	}
	sub_vec->set_uninitialized();  // Nothing to deallocate!
}

template<class Scalar>
void MPIVectorBase<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::MutableSubVectorT<Scalar>* sub_vec )
{
	const Range1D rng = validateRange(rng_in);
	if( rng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		Vector<Scalar>::getSubVector(rng_in,sub_vec);
		return;
	}
	// rng consists of all local data so get it!
	Scalar *localValues = NULL;
	ptrdiff_t stride = 0;
	this->getLocalData(&localValues,&stride);
	sub_vec->initialize(
		rng.lbound()-1                             // globalOffset
		,rng.size()                                // subDim
		,localValues+(rng.lbound()-localOffset_-1) // values
		,stride                                    // stride
		);
}

template<class Scalar>
void MPIVectorBase<Scalar>::commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec )
{
	if( sub_vec->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim() ) {
		// Let the default implementation handle it!
		Vector<Scalar>::commitSubVector(sub_vec);
		return;
	}
	sub_vec->set_uninitialized();  // Nothing to deallocate!
}

// private

template<class Scalar>
Range1D MPIVectorBase<Scalar>::validateRange( const Range1D &rng_in ) const
{
	update_cache();
	const Range1D rng = RangePack::full_range(rng_in,1,globalDim_);
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.lbound() < 1 || globalDim_ < rng.ubound(), std::invalid_argument
		,"EpetraVector::getLocalData(...): Error, the range ["<<rng.lbound()<<","<<rng.ubound()<<"] is not "
		"in the range [1,"<<globalDim_<<"]!"
		);
#endif
	return rng;
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_BASE_HPP
