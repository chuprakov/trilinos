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
// LinearOpWithSolveIterDecl.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_DECL_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_DECL_HPP

#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace TSFCore {
namespace Nonlin {

///
template<class Scalar>
struct LinearOpWithSolveIterState {
	Teuchos::RefCountPtr<const LinearOp<Scalar> >                            M;
	ETransp                                                                       M_trans;
	Teuchos::RefCountPtr<const Solvers::IterativeLinearSolver<Scalar> >      solver;
	Teuchos::RefCountPtr<Solvers::ConvergenceTester<Scalar> >                convTester;
	Teuchos::RefCountPtr<const LinearOp<Scalar> >                            M_tilde_left_inv;
	ETransp                                                                       M_tilde_left_inv_trans;
	Teuchos::RefCountPtr<const LinearOp<Scalar> >                            M_tilde_right_inv;
	ETransp                                                                       M_tilde_right_inv_trans;
};

///
/** Implementation of <tt>LinearOpWithSolve</tt> using an iterative linear solver.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class LinearOpWithSolveIter : virtual public LinearOpWithSolve<Scalar> {
public:
	
	/** @name Constructors / initializers / accessors */
	//@{

	/// Stream that trace information will be sent to
	STANDARD_NONCONST_COMPOSITION_MEMBERS( std::ostream, trace_out )

	///
	LinearOpWithSolveIter();

	///
	LinearOpWithSolveIter(
		const Teuchos::RefCountPtr<const LinearOp<Scalar> >                          &M
		,ETransp                                                                          M_trans
		,const Teuchos::RefCountPtr<const Solvers::IterativeLinearSolver<Scalar> >   &solver
		,const Teuchos::RefCountPtr<Solvers::ConvergenceTester<Scalar> >             &convTester             = Teuchos::null
		,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_left_inv       = Teuchos::null
		,ETransp                                                                          M_tilde_left_inv_trans  = NOTRANS
		,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_right_inv      = Teuchos::null
		,ETransp                                                                          M_tilde_right_inv_trans = NOTRANS
		);

	///
	void initialize(
		const Teuchos::RefCountPtr<const LinearOp<Scalar> >                          &M
		,ETransp                                                                          M_trans
		,const Teuchos::RefCountPtr<const Solvers::IterativeLinearSolver<Scalar> >   &solver
		,const Teuchos::RefCountPtr<Solvers::ConvergenceTester<Scalar> >             &convTester             = Teuchos::null
		,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_left_inv       = Teuchos::null
		,ETransp                                                                          M_tilde_left_inv_trans  = NOTRANS
		,const Teuchos::RefCountPtr<const LinearOp<Scalar> >                         &M_tilde_right_inv      = Teuchos::null
		,ETransp                                                                          M_tilde_right_inv_trans = NOTRANS
		);
	
	///
	LinearOpWithSolveIterState<Scalar> setUninitialized();

	//@}

	/** @name Overridden from OpBase */
	//@{
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > domain() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > range() const;
	///
	bool opSupported(ETransp M_trans) const;
	//@}

	/** @name Overridden from LinearOp */
	//@{
	///
	void apply(
		const ETransp             M_trans
		,const Vector<Scalar>     &x
		,Vector<Scalar>           *y
		,const Scalar             alpha
		,const Scalar             beta
		) const;
	///
	void apply(
		const ETransp               M_trans
		,const MultiVector<Scalar>  &X
		,MultiVector<Scalar>        *Y
		,const Scalar               alpha
		,const Scalar               beta
		) const;
	//@}

	/** @name Overridden from LinearOpWithSolve */
	//@{
	///
	void solve(
		const ETransp                          M_trans
		,const Vector<Scalar>                  &y
		,Vector<Scalar>                        *x
		,Solvers::ConvergenceTester<Scalar>    *convTester
		) const;
	///
	void solve(
		const ETransp                          M_trans
		,const MultiVector<Scalar>             &Y
		,MultiVector<Scalar>                   *X
		,const Scalar                          alpha
		,Solvers::ConvergenceTester<Scalar>    *convTester
		) const;
	///
	Teuchos::RefCountPtr<const LinearOpWithSolve<Scalar> > clone_lows() const;
	///
	Teuchos::RefCountPtr<const LinearOp<Scalar> > preconditioner() const;
	//@}

private:

#ifdef DOXYGEN_COMPILE
	LinearOp<Scalar>                              *M;
	Solvers::IterativeLinearSolver<Scalar>        *solver;
	Solvers::ConvergenceTester<Scalar>            *convTester;
	LinearOp<<Scalar>                             *M_tilde_left_inv;
	LinearOp<<Scalar>                             *M_tilde_right_inv;
#else
	LinearOpWithSolveIterState<Scalar>            state_;
#endif

}; // class LinearOpWithSolveIter

} // namespace Nonlin
} // namespace TSFCore

#endif	// TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_ITER_DECL_HPP
