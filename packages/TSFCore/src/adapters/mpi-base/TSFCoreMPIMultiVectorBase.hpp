// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIMultiVectorBase.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_HPP

#include <vector>

#include "TSFCoreMPIMultiVectorBaseDecl.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

template<class Scalar>
MPIMultiVectorBase<Scalar>::MPIMultiVectorBase()
	:in_applyOp_(false)
{}

// Overridden form OpBase

template<class Scalar>
MemMngPack::ref_count_ptr< const VectorSpace<Scalar> >
MPIMultiVectorBase<Scalar>::range() const
{
	return mpiSpace();
}

// Overridden from MultiVector

template<class Scalar>
void MPIMultiVectorBase<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &primary_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_objs[]
	,const Index                    primary_first_ele
	,const Index                    primary_sub_dim
	,const Index                    primary_global_offset
	,const Index                    secondary_first_ele
	,const Index                    secondary_sub_dim
	) const
{
	// ToDo: Provide a specialized implementation!
	MultiVector<Scalar>::applyOp(
		primary_op
		,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
		,reduct_objs
		,primary_first_ele,primary_sub_dim,primary_global_offset
		,secondary_first_ele,secondary_sub_dim
		);
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_HPP
