// //////////////////////////////////////////////////////////////////////////////////
// LinearOpWithSolveAztecOO.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_AZTECOO_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_AZTECOO_HPP

#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "StandardCompositionMacros.hpp"
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

  ///
  typedef double Scalar;
	
	/** @name Constructors / initializers / accessors */
	//@{

	/// Stream that trace to which information will be sent
	STANDARD_NONCONST_COMPOSITION_MEMBERS( std::ostream, trace_out );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxIter );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, relTol );

	///
	/** Construct uninitialized but with default option values.
	 *
	 * Note, these defaults where taken from
	 * NOX::Epetra::Group::applyJacobianInverse(...) on 2004/01/19.
	 */
 	LinearOpWithSolveAztecOO(
	 	const int      maxIter = 400
		,const double  relTol  = 1e-6
		);

	///
	LinearOpWithSolveAztecOO(
		const Teuchos::RefCountPtr<Epetra_Operator>                          &Op
		,ETransp                                                             Op_trans
		,const Teuchos::RefCountPtr<AztecOO>                                 &solver
		,const Teuchos::RefCountPtr<Epetra_Operator>                         &Prec               = Teuchos::null
		,ETransp                                                             Prec_trans          = NOTRANS
    ,bool                                                                adjointSupported    = false
		);

	///
	void initialize(
		const Teuchos::RefCountPtr<Epetra_Operator>                          &Op
		,ETransp                                                             Op_trans
		,const Teuchos::RefCountPtr<AztecOO>                                 &solver
		,const Teuchos::RefCountPtr<Epetra_Operator>                         &Prec               = Teuchos::null
		,ETransp                                                             Prec_trans          = NOTRANS
    ,bool                                                                adjointSupported    = false
		);
	
	///
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_Operator>                                *Op                  = NULL
		,ETransp                                                             *Op_trans            = NULL
		,Teuchos::RefCountPtr<AztecOO>                                       *solver              = NULL
		,Teuchos::RefCountPtr<Epetra_Operator>                               *Prec                = NULL
		,ETransp                                                             *Prec_trans          = NULL
    ,bool                                                                *adjointSupported    = NULL
    );

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
  Teuchos::RefCountPtr<Epetra_Operator>         Op_;
	ETransp                                       Op_trans_;
  EpetraLinearOp                                tsfcore_Op_;
  Teuchos::RefCountPtr<AztecOO>                 solver_;
  Teuchos::RefCountPtr<Epetra_Operator>         Prec_;
  ETransp                                       Prec_trans_;
  EpetraLinearOp                                tsfcore_Prec_;
  bool                                          adjointSupported_;
#endif

}; // class LinearOpWithSolveAztecOO

} // namespace Nonlin
} // namespace TSFCore

#endif	// TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_AZTECOO_HPP
