// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearSolveOp.hpp

#ifndef TSFCORE_NONLIN_LINEAR_SOLVE_OP_HPP
#define TSFCORE_NONLIN_LINEAR_SOLVE_OP_HPP

#include <assert.h>

#include "TSFCoreNonlinLinearSolveOpDecl.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVectorStdOps.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
ENonsingStatus LinearSolveOp<Scalar>::nonsingStatus() const
{
	return OP_SINGULARITY_UNKNOWN;
}

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearSolveOp<Scalar> >
LinearSolveOp<Scalar>::clone_lso() const
{
	return MemMngPack::null;
}

template<class Scalar>
void LinearSolveOp<Scalar>::solve(
	const ETransp                          M_trans
	,const MultiVector<Scalar>             &Y
	,MultiVector<Scalar>                   *X
	,const Scalar                          alpha
	,Solvers::ConvergenceTester<Scalar>    *convTester
	) const
{
	namespace mmp = MemMngPack;
	Mp_MtM_assert_compatibility(X,NOTRANS,*this,M_trans==NOTRANS?TRANS:NOTRANS,Y,NOTRANS);
	const VectorSpace<Scalar> &space_mv_rows = *Y.domain();
	const Index               num_mv_cols    = space_mv_rows.dim();
	//
	// Here we will solve the linear systems one at a time as:
	//
	//    op(M)*X(j) = alpha*Y(j))
	//
	bool scale_y = (alpha != 1.0);
	mmp::ref_count_ptr<Vector<Scalar> > tmp = ( scale_y ? Y.range()->createMember() : mmp::null );
	for( Index j = 1; j <= num_mv_cols; ++j ) {
		// get tmp = alpha*Y.col(j) (but only scale if alpha != 1.0)
		if( scale_y ) {
			*tmp = *Y.col(j);
			Vt_S(tmp.get(),alpha);
		}
		else {
			tmp = mmp::rcp_const_cast<Vector<Scalar> >(Y.col(j));
		}
		// Solve op(M)*X(j) = alpha*Y(j)
		this->solve(M_trans,*tmp,X->col(j).get(),convTester);
	}
}

} // namespace Nonlin
} // namespace TSFCore

#endif /// TSFCORE_NONLIN_LINEAR_SOLVE_OP_HPP
