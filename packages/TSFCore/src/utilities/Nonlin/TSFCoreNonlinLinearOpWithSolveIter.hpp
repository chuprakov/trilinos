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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearOpWithSolveIter.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_HPP

#include "TSFCoreNonlinLinearOpWithSolveIterDecl.hpp"
#include "TSFCoreSolversIterativeLinearSolver.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Nonlin {

// Constructors / initializers / accessors

template<class Scalar>
LinearOpWithSolveIter<Scalar>::LinearOpWithSolveIter()
{}

template<class Scalar>
LinearOpWithSolveIter<Scalar>::LinearOpWithSolveIter(
	const Teuchos::RefCountPtr<const LinearOp<Scalar> >                          &M
	,ETransp                                                                     M_trans
	,const Teuchos::RefCountPtr<const Solvers::IterativeLinearSolver<Scalar> >   &solver
	,const Teuchos::RefCountPtr<Solvers::ConvergenceTester<Scalar> >             &convTester
	,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_left_inv
	,ETransp                                                                     M_tilde_left_inv_trans
	,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_right_inv
	,ETransp                                                                     M_tilde_right_inv_trans
	)
{
	initialize(M,M_trans,solver,convTester,M_tilde_left_inv,M_tilde_left_inv_trans,M_tilde_right_inv,M_tilde_right_inv_trans);
}

template<class Scalar>
void LinearOpWithSolveIter<Scalar>::initialize(
	const Teuchos::RefCountPtr<const LinearOp<Scalar> >                          &M
	,ETransp                                                                     M_trans
	,const Teuchos::RefCountPtr<const Solvers::IterativeLinearSolver<Scalar> >   &solver
	,const Teuchos::RefCountPtr<Solvers::ConvergenceTester<Scalar> >             &convTester
	,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_left_inv
	,ETransp                                                                     M_tilde_left_inv_trans
	,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_right_inv
	,ETransp                                                                     M_tilde_right_inv_trans
	)
{
#ifdef _DEBUG
	const char func_name[] = "LinearOpWithSolveIter<Scalar>::initialize(...)";
	TEST_FOR_EXCEPTION(M.get()==NULL,std::invalid_argument,func_name<<": Error!");
	TEST_FOR_EXCEPTION(solver.get()==NULL,std::invalid_argument,func_name<<": Error!");
	if(solver->adjointRequired()) {
		const bool adjoint_supported = M->opSupported(not_trans(M_trans));
		TEST_FOR_EXCEPTION(
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
		TEST_FOR_EXCEPTION(
			!(domain_compatible && range_compatible), Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_left_inv) do are not compatible!");
		if(solver->adjointRequired()) {
			const bool adjoint_supported = M_tilde_left_inv->opSupported(not_trans(M_tilde_left_inv_trans));
			TEST_FOR_EXCEPTION(
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
		TEST_FOR_EXCEPTION(
			!(domain_compatible && range_compatible), Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_right_inv) do are not compatible!");
		if(solver->adjointRequired()) {
			const bool adjoint_supported = M_tilde_right_inv->opSupported(not_trans(M_tilde_right_inv_trans));
			TEST_FOR_EXCEPTION(
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

	state_.M                       = Teuchos::null;
	state_.M_trans                 = NOTRANS;
	state_.solver                  = Teuchos::null;
	state_.convTester              = Teuchos::null;
	state_.M_tilde_left_inv        = Teuchos::null;
	state_.M_tilde_left_inv_trans  = NOTRANS;
	state_.M_tilde_right_inv       = Teuchos::null;
	state_.M_tilde_right_inv_trans = NOTRANS;

	return state_tmp;
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinearOpWithSolveIter<Scalar>::domain() const
{
	return (state_.M_trans == NOTRANS ? state_.M->domain() : state_.M->range() );
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinearOpWithSolveIter<Scalar>::range() const
{
	return (state_.M_trans == NOTRANS ? state_.M->range() : state_.M->domain() );
}

template<class Scalar>
bool LinearOpWithSolveIter<Scalar>::opSupported(ETransp M_trans) const
{
	return state_.M->opSupported(M_trans);
}

// Overridden from LinearOp

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
#ifdef TSFCORE_VECTOR_DERIVE_FROM_MULTI_VECTOR
	solve(M_trans,static_cast<const MultiVector<Scalar>&>(y),static_cast<MultiVector<Scalar>*>(x),1.0,convTester);
#else
	const MultiVectorCols<Scalar>  Y(Teuchos::rcp(const_cast<Vector<Scalar>*>(&y),false));
	MultiVectorCols<Scalar>        X(Teuchos::rcp(x,false));
	solve(M_trans,Y,&X,1.0,convTester);
#endif
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
	if(convTester) convTester->attach(state_.convTester);
	Solvers::SolveReturn
		solve_return = state_.solver->solve(
			*state_.M, trans_trans( M_trans, state_.M_trans )
			,Y, X, alpha
			,Solvers::DEFAULT_MAX_ITER
			,convTester ? convTester : state_.convTester.get()
			,state_.M_tilde_left_inv.get(),  trans_trans( M_trans, state_.M_tilde_left_inv_trans )
			,state_.M_tilde_right_inv.get(), trans_trans( M_trans, state_.M_tilde_right_inv_trans )
			);
	if(convTester) convTester->attach(Teuchos::null);
	switch(solve_return.solve_status) {
		case Solvers::SOLVED_TO_TOL: {
			if(get_trace_out().get())
				trace_out() << "\nLinear system(s) solved to tolerance in num_iter = "<<solve_return.num_iter<<" iterations!\n";
			// Great! we solved it!
			break;
		}
		case Solvers::MAX_ITER_EXCEEDED: {
			TEST_FOR_EXCEPTION(
				true, Solvers::Exceptions::FailureToConverge
				,"LinearOpWithSolveIter<Scalar>::solve(...): Error, num_iter = "<<solve_return.num_iter<<" iterations were "
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
Teuchos::RefCountPtr<const LinearOpWithSolve<Scalar> >
LinearOpWithSolveIter<Scalar>::clone_lows() const
{
	namespace mmp = MemMngPack;
	if(state_.M.get()) {
		return Teuchos::rcp(
			new LinearOpWithSolveIter(
				state_.M->clone()
				,state_.M_trans
				,state_.solver->clone()
				,(state_.convTester.get() ? state_.convTester->clone() : Teuchos::null )
				,(state_.M_tilde_left_inv.get() ? state_.M_tilde_left_inv->clone() : Teuchos::null)
				,state_.M_tilde_left_inv_trans
				,(state_.M_tilde_right_inv.get() ? state_.M_tilde_right_inv->clone() : Teuchos::null)
				,state_.M_tilde_right_inv_trans
				) );
	}
	return Teuchos::rcp( new LinearOpWithSolveIter() );
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
LinearOpWithSolveIter<Scalar>::preconditioner() const
{
	TEST_FOR_EXCEPT(true); // ToDo: Return a composite operator for M_tilde_left_inv*M_tilde_right_inv
	return state_.M_tilde_left_inv;
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_HPP
