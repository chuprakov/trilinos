// ////////////////////////////////////////////////
// TSFCoreExplInst.cpp

#include "TSFCoreVectorSpaceFactory.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreOpBase.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreSerialVectorSpaceFactory.hpp"
#include "TSFCoreSerialVectorSpaceBase.hpp"
#include "TSFCoreSerialVectorSpace.hpp"
#include "TSFCoreSerialVectorBase.hpp"
#include "TSFCoreSerialVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"

/* 
#include "Solvers/Norm.hpp"
#include "Solvers/SolverState.hpp"
#include "Solvers/ConvergenceTester.hpp"
#include "Solvers/IterativeLinearSolver.hpp"
#include "Nonlin/LinearSolveOp.hpp"
#include "Nonlin/LinearOpWithSolve.hpp"
#include "Nonlin/NonlinearProblem.hpp"
#include "Nonlin/NonlinearProblemFirstOrder.hpp"
*/

namespace TSFCore {

template VectorSpaceFactory<RTOp_value_type>;
template VectorSpace<RTOp_value_type>;
template Vector<RTOp_value_type>;
template OpBase<RTOp_value_type>;
template LinearOp<RTOp_value_type>;
template MultiVector<RTOp_value_type>;
template SerialVectorSpaceFactory<RTOp_value_type>;
template SerialVectorSpaceBase<RTOp_value_type>;
template SerialVectorSpace<RTOp_value_type>;
template SerialVectorBase<RTOp_value_type>;
template SerialVector<RTOp_value_type>;
template MultiVectorCols<RTOp_value_type>;

template void applyOp<RTOp_value_type>(
	const RTOpPack::RTOpT<RTOp_value_type>    &op
	,const size_t                             num_vecs
	,const Vector<RTOp_value_type>*           vecs[]
	,const size_t                             num_targ_vecs
	,Vector<RTOp_value_type>*                 targ_vecs[]
	,RTOp_ReductTarget                        reduct_obj
	,const Index                              first_ele
	,const Index                              sub_dim
	,const Index                              global_offset
	);

  /*
namespace Solvers {

template Norm<RTOp_value_type>;
template ConvergenceTester<RTOp_value_type>;
template IterativeLinearSolver<RTOp_value_type>;

} // namespace Solvers

namespace Nonlin {

template LinearSolveOp<RTOp_value_type>;
template LinearOpWithSolve<RTOp_value_type>;
template NonlinearProblem<RTOp_value_type>;
template NonlinearProblemFirstOrder<RTOp_value_type>;

} // namespace Nonlin

  */
}  // namespace TSFCore
