// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreLinearOp.hpp

#ifndef TSFCORE_LINEAR_OP_HPP
#define TSFCORE_LINEAR_OP_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreLinearOpDecl.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreAssertOp.hpp"

namespace TSFCore {

// Virtual functions with default implemenations

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearOp<Scalar> > 
LinearOp<Scalar>::clone() const
{
	return MemMngPack::null;
}

template<class Scalar>
void LinearOp<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	const VectorSpace<Scalar> &space_mv_rows = *Y->domain();
	const Index               num_mv_cols    = space_mv_rows.dim();
	for( Index j = 1; j <= num_mv_cols; ++j )
		this->apply(M_trans,*X.col(j),Y->col(j).get(),alpha,beta);
}

}	// end namespace TSFCore

#endif // TSFCORE_LINEAR_OP_HPP
