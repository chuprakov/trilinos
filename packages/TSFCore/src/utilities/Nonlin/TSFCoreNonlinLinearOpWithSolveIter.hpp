// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearOpWithSolveIter.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_HPP

#include "TSFCoreNonlinLinearOpWithSolveIterDecl.hpp"
#include "TSFCoreSolversIterativeLinearSolver.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "ThrowException.hpp"

namespace TSFCore {
namespace Nonlin {

// Constructors / initializers / accessors

template<class Scalar>
LinearOpWithSolveIter<Scalar>::LinearOpWithSolveIter()
{}

template<class Scalar>
LinearOpWithSolveIter<Scalar>::LinearOpWithSolveIter(
	const MemMngPack::ref_count_ptr<const LinearOp<Scalar> >                          &M
	,ETransp                                                                          M_trans
	,const MemMngPack::ref_count_ptr<const Solvers::IterativeLinearSolver<Scalar> >   &solver
	,const MemMngPack::ref_count_ptr<Solvers::ConvergenceTester<Scalar> >             &convTester
	,const MemMngPack::ref_count_ptr<const LinearOp<Scalar> >                         &M_tilde_left_inv
	,ETransp                                                                          M_tilde_left_inv_trans
	,const MemMngPack::ref_count_ptr<const LinearOp<Scalar> >                         &M_tilde_right_inv
	,ETransp                                                                          M_tilde_right_inv_trans
	)
{
	initialize(M,M_trans,solver,convTester,M_tilde_left_inv,M_tilde_left_inv_trans,M_tilde_right_inv,M_tilde_right_inv_trans);
}

template<class Scalar>
void LinearOpWithSolveIter<Scalar>::initialize(
	const MemMngPack::ref_count_ptr<const LinearOp<Scalar> >                          &M
	,ETransp                                                                          M_trans
	,const MemMngPack::ref_count_ptr<const Solvers::IterativeLinearSolver<Scalar> >   &solver
	,const MemMngPack::ref_count_ptr<Solvers::ConvergenceTester<Scalar> >             &convTester
	,const MemMngPack::ref_count_ptr<const LinearOp<Scalar> >                         &M_tilde_left_inv
	,ETransp                                                                          M_tilde_left_inv_trans
	,const MemMngPack::ref_count_ptr<const LinearOp<Scalar> >                         &M_tilde_right_inv
	,ETransp                                                                          M_tilde_right_inv_trans
	)
{
#ifdef _DEBUG
	const char func_name[] = "LinearOpWithSolveIter<Scalar>::initialize(...)";
	THROW_EXCEPTION(M.get()==NULL,std::invalid_argument,func_name<<": Error!");
	THROW_EXCEPTION(solver.get()==NULL,std::invalid_argument,func_name<<": Error!");
	if(solver->adjointRequired()) {
		const bool adjoint_supported = M->opSupported(not_trans(M_trans));
		THROW_EXCEPTION(
			!adjoint_supported,std::invalid_argument
			,func_name<<": Error, solver requires support for both non-transposed and transposed ops and M does not comply!");
	}
	const VectorSpace<Scalar> &opM_domain = ( M_trans==NOTRANS ? *M->domain() : *M->range()  );
	const VectorSpace<Scalar> &opM_range  = ( M_trans==NOTRANS ? *M->range()  : *M->domain() );
	if(M_tilde_left_inv.get()) {
		const VectorSpace<Scalar> &opM_tilde_left_inv_domain = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->domain() : *M_tilde_left_inv->range()  );
		const VectorSpace<Scalar> &opM_tilde_left_inv_range  = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->range()  : *M_tilde_left_inv->domain() );
		const bool
			domain_compatible = opM_domain.isCompatible(opM_tilde_left_inv_range),
			range_compatible =  opM_range.isCompatible(opM_tilde_left_inv_domain);
		THROW_EXCEPTION(
			!(domain_compatible && range_compatible), Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_left_inv) do are not compatible!");
		if(solver->adjointRequired()) {
			const bool adjoint_supported = M_tilde_left_inv->opSupported(not_trans(M_tilde_left_inv_trans));
			THROW_EXCEPTION(
				!adjoint_supported,std::invalid_argument
				,func_name<<": Error, solver requires support for both non-transposed and transposed ops and M_tilde_left_inv does not comply!");
		}
	}
	if(M_tilde_right_inv.get()) {
		const VectorSpace<Scalar> &opM_tilde_right_inv_domain = ( M_tilde_right_inv_trans==NOTRANS ? *M_tilde_right_inv->domain() : *M_tilde_right_inv->range()  );
		const VectorSpace<Scalar> &opM_tilde_right_inv_range  = ( M_tilde_right_inv_trans==NOTRANS ? *M_tilde_right_inv->range()  : *M_tilde_right_inv->domain() );
		const bool
			domain_compatible = opM_domain.isCompatible(opM_tilde_right_inv_range),
			range_compatible =  opM_range.isCompatible(opM_tilde_right_inv_domain);
		THROW_EXCEPTION(
			!(domain_compatible && range_compatible), Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_right_inv) do are not compatible!");
		if(solver->adjointRequired()) {
			const bool adjoint_supported = M_tilde_right_inv->opSupported(not_trans(M_tilde_right_inv_trans));
			THROW_EXCEPTION(
				!adjoint_supported,std::invalid_argument
				,func_name<<": Error, solver requires support for both non-transposed and transposed ops and M_tilde_right_inv does not comply!");
		}
	}
#endif
	state_.M                       = M;
	state_.M_trans                 = M_trans;
	state_.solver                  = solver;
	state_.convTester              = convTester;
	state_.M_tilde_left_inv        = M_tilde_left_inv;
	state_.M_tilde_left_inv_trans  = M_tilde_left_inv_trans;
	state_.M_tilde_right_inv       = M_tilde_right_inv;
	state_.M_tilde_right_inv_trans = M_tilde_right_inv_trans;
}
	
template<class Scalar>
LinearOpWithSolveIterState<Scalar>
LinearOpWithSolveIter<Scalar>::setUninitialized()
{
	namespace mmp = MemMngPack;

	LinearOpWithSolveIterState<Scalar> state_tmp = state_;

	state_.M                       = mmp::null;
	state_.M_trans                 = NOTRANS;
	state_.solver                  = mmp::null;
	state_.convTester              = mmp::null;
	state_.M_tilde_left_inv        = mmp::null;
	state_.M_tilde_left_inv_trans  = NOTRANS;
	state_.M_tilde_right_inv       = mmp::null;
	state_.M_tilde_right_inv_trans = NOTRANS;

	return state_tmp;
}

// Overridden from LinearOp

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >
LinearOpWithSolveIter<Scalar>::domain() const
{
	return (state_.M_trans == NOTRANS ? state_.M->domain() : state_.M->range() );
}

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >
LinearOpWithSolveIter<Scalar>::range() const
{
	return (state_.M_trans == NOTRANS ? state_.M->range() : state_.M->domain() );
}

template<class Scalar>
bool LinearOpWithSolveIter<Scalar>::opSupported(ETransp M_trans) const
{
	return state_.M->opSupported(M_trans);
}

template<class Scalar>
void LinearOpWithSolveIter<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	state_.M->apply(trans_trans(M_trans,state_.M_trans),x,y,alpha,beta);
}

template<class Scalar>
void LinearOpWithSolveIter<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	state_.M->apply(trans_trans(M_trans,state_.M_trans),X,Y,alpha,beta);
}

// Overridden from LinearOpWithSolve

template<class Scalar>
void LinearOpWithSolveIter<Scalar>::solve(
	const ETransp                        M_trans
	,const Vector<Scalar>                &y
	,Vector<Scalar>                      *x
	,Solvers::ConvergenceTester<Scalar>  *convTester
	) const
{
	namespace mmp = MemMngPack;
	const MultiVectorCols<Scalar>  Y(mmp::rcp(const_cast<Vector<Scalar>*>(&y),false));
	MultiVectorCols<Scalar>        X(mmp::rcp(x,false));
	solve(M_trans,Y,&X,1.0,convTester);
}

template<class Scalar>
void LinearOpWithSolveIter<Scalar>::solve(
	const ETransp                         M_trans
	,const MultiVector<Scalar>            &Y
	,MultiVector<Scalar>                  *X
	,const Scalar                         alpha
	,Solvers::ConvergenceTester<Scalar>   *convTester
	) const
{
	if(get_trace_out().get())
		trace_out()
			<< "\n*** Entering LinearOpWithSolveIter<Scalar>::solve(...):...\n"
			<< "\nUsing a linear solver of type \'" << typeid(*state_.solver).name() << "\' ...\n";
	Solvers::SolveReturn
		solve_return = state_.solver->solve(
			*state_.M, trans_trans( M_trans, state_.M_trans )
			,Y, X, alpha
			,Solvers::DEFAULT_MAX_ITER
			,convTester ? convTester : state_.convTester.get()
			,state_.M_tilde_left_inv.get(),  trans_trans( M_trans, state_.M_tilde_left_inv_trans )
			,state_.M_tilde_right_inv.get(), trans_trans( M_trans, state_.M_tilde_right_inv_trans )
			);
	switch(solve_return.solve_status) {
		case Solvers::SOLVED_TO_TOL: {
			if(get_trace_out().get())
				trace_out() << "\nLinear system(s) solved to tolerance in num_iter = "<<solve_return.num_iter<<" iterations!\n";
			// Great! we solved it!
			break;
		}
		case Solvers::MAX_ITER_EXCEEDED: {
			THROW_EXCEPTION(
				true, Solvers::Exceptions::FailureToConverge
				,"LinearOpWithSolveIter<Scalar>::solve(...): Error, num_iter = "<<solve_return.num_iter<<" iterations where "
				"performed by the solver object and exceeded the maximum number!"
				);
			break;
		}
		default: {
			assert(0);
		}
	}
	if(get_trace_out().get())
		trace_out() << "\n*** Leaving LinearOpWithSolveIter<Scalar>::solve(...) ...\n\n";
}

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearOpWithSolve<Scalar> >
LinearOpWithSolveIter<Scalar>::clone_lows() const
{
	namespace mmp = MemMngPack;
	if(state_.M.get()) {
		return mmp::rcp(
			new LinearOpWithSolveIter(
				state_.M->clone()
				,state_.M_trans
				,state_.solver->clone()
				,(state_.convTester.get() ? state_.convTester->clone() : mmp::null )
				,(state_.M_tilde_left_inv.get() ? state_.M_tilde_left_inv->clone() : mmp::null)
				,state_.M_tilde_left_inv_trans
				,(state_.M_tilde_right_inv.get() ? state_.M_tilde_right_inv->clone() : mmp::null)
				,state_.M_tilde_right_inv_trans
				) );
	}
	return mmp::rcp( new LinearOpWithSolveIter() );
}

template<class Scalar>
MemMngPack::ref_count_ptr<const LinearOp<Scalar> >
LinearOpWithSolveIter<Scalar>::preconditioner() const
{
	assert(0); // ToDo: Return a composite operator for M_tilde_left_inv*M_tilde_right_inv
	return state_.M_tilde_left_inv;
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_HPP
