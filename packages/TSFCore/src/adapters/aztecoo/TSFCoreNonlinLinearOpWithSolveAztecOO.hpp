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
// LinearOpWithSolveAztecOO.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_AZTECOO_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_AZTECOO_HPP

#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "Epetra_LinearSystemScaler.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "AztecOO.h"

namespace TSFCore {
namespace Nonlin {

///
/** Implementation of <tt>LinearOpWithSolve</tt> using AztecOO and Epetra.
 *
 * ToDo: Finish documentation!
 */
class LinearOpWithSolveAztecOO : virtual public LinearOpWithSolve<double> {
public:

	/** @name Public types */
	//@{

  ///
  typedef double Scalar;

	//@}
	
	/** @name Constructors / initializers / accessors */
	//@{

	/// Stream that trace to which information will be sent
	STANDARD_NONCONST_COMPOSITION_MEMBERS( std::ostream, trace_out )
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxIter )
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, relTol )
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, minRelTol )

	///
	/** Give non-const access to the linear system scaler object so that
	 * clients can set options.
	 */
	Epetra::LinearSystemScaler& linearSystemScaler();

	///
	const Epetra::LinearSystemScaler& linearSystemScaler() const;

	///
	/** Construct uninitialized but with default option values.
	 *
	 * Note, these defaults where taken from
	 * NOX::Epetra::Group::applyJacobianInverse(...) on 2004/01/19.
	 */
 	LinearOpWithSolveAztecOO(
	 	const int                    maxIter            = 400
		,const double                relTol             = 1e-6
		,const double                minRelTol          = 1e-2
		);

	///
	~LinearOpWithSolveAztecOO();

	///
	/** Sets up this object with an operator and an optional preconditioner.
	 *
	 * Note: All options must be set before calling this function.
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_Operator>      &Op
		,const ETransp                                   Op_trans
		,const Teuchos::RefCountPtr<AztecOO>             &solver
		,const Teuchos::RefCountPtr<Epetra_Operator>     &Prec               = Teuchos::null
		,const ETransp                                   Prec_trans          = NOTRANS
		,const Epetra::ProductOperator::EApplyMode       Prec_inverse        = Epetra::ProductOperator::APPLY_MODE_APPLY_INVERSE
    ,const bool                                      adjointSupported    = false
		);
	
	///
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_Operator>            *Op                  = NULL
		,ETransp                                         *Op_trans            = NULL
		,Teuchos::RefCountPtr<AztecOO>                   *solver              = NULL
		,Teuchos::RefCountPtr<Epetra_Operator>           *Prec                = NULL
		,ETransp                                         *Prec_trans          = NULL
		,Epetra::ProductOperator::EApplyMode             *Prec_inverse        = NULL
    ,bool                                            *adjointSupported    = NULL
    );

	/// Returns <tt>return.get()!=NULL</tt> if <tt>*this</tt> is initialized!.
	Teuchos::RefCountPtr<Epetra_Operator> Op() const;

	///
	void resetCounters();

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

  Epetra_Operator          *Op;
  AztecOO                  *solver;
  Epetra_Operator          *Prec;

#else

	Epetra::LinearSystemScaler                    linearSystemScaler_;

  Teuchos::RefCountPtr<Epetra_Operator>         Op_;
	ETransp                                       Op_trans_;
  EpetraLinearOp                                tsfcore_Op_;
  Teuchos::RefCountPtr<AztecOO>                 solver_;
  Teuchos::RefCountPtr<Epetra_Operator>         Prec_;
  ETransp                                       Prec_trans_;
	Epetra::ProductOperator::EApplyMode           Prec_inverse_;
  EpetraLinearOp                                tsfcore_Prec_;
	
  bool                                          adjointSupported_;

  Teuchos::RefCountPtr<Epetra_Operator>         fwd_Op_;
  Teuchos::RefCountPtr<Epetra_Operator>         fwd_Prec_;
  Teuchos::RefCountPtr<Epetra_Operator>         adj_Op_;
  Teuchos::RefCountPtr<Epetra_Operator>         adj_Prec_;

	mutable int                                   numFwdSolves_;
	mutable int                                   numFwdLinearIters_;
	mutable double                                fwdLinearCPU_; // seconds
	mutable int                                   numAdjSolves_;
	mutable int                                   numAdjLinearIters_;
	mutable double                                adjLinearCPU_; // seconds

#endif

}; // class LinearOpWithSolveAztecOO

// ////////////////////////////
// Inline members

inline
Epetra::LinearSystemScaler&
LinearOpWithSolveAztecOO::linearSystemScaler()
{
	return linearSystemScaler_;
}

inline
const Epetra::LinearSystemScaler&
LinearOpWithSolveAztecOO::linearSystemScaler() const
{
	return linearSystemScaler_;
}

inline
Teuchos::RefCountPtr<Epetra_Operator>
LinearOpWithSolveAztecOO::Op() const
{
	return Op_;
}

} // namespace Nonlin
} // namespace TSFCore

#endif	// TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_AZTECOO_HPP
