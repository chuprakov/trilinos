// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearOpWithSolveAztecOO.cpp

#include "TSFCoreNonlinLinearOpWithSolveAztecOO.hpp"
#include "TSFCore_get_Epetra_MultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Nonlin {

// Constructors / initializers / accessors

LinearOpWithSolveAztecOO::LinearOpWithSolveAztecOO(
	const int      maxIter
	,const double  relTol
	,const double  minRelTol
	)
	:maxIter_(maxIter)
	,relTol_(relTol)
	,minRelTol_(minRelTol)
{}

LinearOpWithSolveAztecOO::LinearOpWithSolveAztecOO(
  const Teuchos::RefCountPtr<Epetra_Operator>                          &Op
  ,ETransp                                                             Op_trans
  ,const Teuchos::RefCountPtr<AztecOO>                                 &solver
  ,const Teuchos::RefCountPtr<Epetra_Operator>                         &Prec
  ,ETransp                                                             Prec_trans
  ,bool                                                                adjointSupported
	)
{
	initialize(Op,Op_trans,solver,Prec,Prec_trans,adjointSupported);
}

void LinearOpWithSolveAztecOO::initialize(
  const Teuchos::RefCountPtr<Epetra_Operator>                          &Op
  ,ETransp                                                             Op_trans
  ,const Teuchos::RefCountPtr<AztecOO>                                 &solver
  ,const Teuchos::RefCountPtr<Epetra_Operator>                         &Prec
  ,ETransp                                                             Prec_trans
  ,bool                                                                adjointSupported
	)
{
#ifdef _DEBUG
	const char func_name[] = "LinearOpWithSolveAztecOO::initialize(...)";
	TEST_FOR_EXCEPTION(Op.get()==NULL,std::invalid_argument,func_name<<": Error!");
	TEST_FOR_EXCEPTION(solver.get()==NULL,std::invalid_argument,func_name<<": Error!");
  // ToDo: Check vector spaces!
#endif
  Op_ = Op;
  Op_trans_ = Op_trans;
  tsfcore_Op_.initialize(Op,Op_trans);
  solver_ = solver;
  if(Prec.get()) tsfcore_Prec_.initialize(Prec,Prec_trans);
  Prec_ = Prec;
  Prec_trans_ = Prec_trans;
  adjointSupported_ = adjointSupported;
}
	
void LinearOpWithSolveAztecOO::setUninitialized(
  Teuchos::RefCountPtr<Epetra_Operator>                                *Op
  ,ETransp                                                             *Op_trans
  ,Teuchos::RefCountPtr<AztecOO>                                       *solver
  ,Teuchos::RefCountPtr<Epetra_Operator>                               *Prec
  ,ETransp                                                             *Prec_trans
  ,bool                                                                *adjointSupported
  )
{
  Op_ = Teuchos::null;
  Op_trans_ = NOTRANS;
  tsfcore_Op_.setUninitialized(Op,Op_trans);
  if(solver) *solver = solver_; solver_ = Teuchos::null;
  Prec_ = Teuchos::null;
  Prec_trans_ = NOTRANS;
  tsfcore_Prec_.setUninitialized(Prec,Prec_trans);
  if(adjointSupported) *adjointSupported = adjointSupported_;
}

// Overridden from OpBase

Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinearOpWithSolveAztecOO::domain() const
{
  return tsfcore_Op_.domain();
}

Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinearOpWithSolveAztecOO::range() const
{
  return tsfcore_Op_.range();
}

bool LinearOpWithSolveAztecOO::opSupported(ETransp M_trans) const
{
  return ( M_trans == NOTRANS ) || ( M_trans == TRANS && adjointSupported_ );
}

// Overridden from LinearOp

void LinearOpWithSolveAztecOO::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
  tsfcore_Op_.apply(M_trans,x,y,alpha,beta);
}

void LinearOpWithSolveAztecOO::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
  tsfcore_Op_.apply(M_trans,X,Y,alpha,beta);
}

// Overridden from LinearOpWithSolve

void LinearOpWithSolveAztecOO::solve(
	const ETransp                        M_trans
	,const Vector<Scalar>                &y
	,Vector<Scalar>                      *x
	,Solvers::ConvergenceTester<Scalar>  *convTester
	) const
{
	const MultiVectorCols<Scalar>  Y(Teuchos::rcp(const_cast<Vector<Scalar>*>(&y),false));
	MultiVectorCols<Scalar>        X(Teuchos::rcp(x,false));
	solve(M_trans,Y,&X,1.0,convTester);
}

void LinearOpWithSolveAztecOO::solve(
	const ETransp                         M_trans
	,const MultiVector<Scalar>            &Y_in
	,MultiVector<Scalar>                  *X_inout
	,const Scalar                         alpha
	,Solvers::ConvergenceTester<Scalar>   *convTester
	) const
{
	if(get_trace_out().get())
		trace_out()
			<< "\n*** Entering LinearOpWithSolveAztecOO::solve(...):...\n";
  //
  // 1/alpha*op(M)*X = Y
	//
	// Get Epetra_MultiVector objects for the arguments
	//
	Teuchos::RefCountPtr<const Epetra_MultiVector>
		Y = get_Epetra_MultiVector(
			M_trans==NOTRANS ? *tsfcore_Op_.epetraRange() : *tsfcore_Op_.epetraDomain()
			,Teuchos::rcp(&Y_in,false)
			);
	Teuchos::RefCountPtr<Epetra_MultiVector>
		X = get_Epetra_MultiVector(
			M_trans==NOTRANS ? *tsfcore_Op_.epetraDomain() : *tsfcore_Op_.epetraRange()
			,Teuchos::rcp(X_inout,false)
			);
  //
  // Create the linear problem
  //
  if(get_trace_out().get())
    trace_out() << "\nUsing an operator of type \'" << typeid(*Op_).name() << "\'\n";
  Op_->SetUseTranspose(trans_trans(Op_trans_,M_trans)==TRANS); // Must set mode just before use!
  Epetra_LinearProblem
    lp(
      &*Op_
      ,&*X
      ,const_cast<Epetra_MultiVector*>(&*Y) // Epetra_LinearProblem should be fixed!
      );
  lp.SetPDL(::hard); // Set a hard problem by default
  //
  // Create new AztecOO solver and copy options from externally given solver
  //
  AztecOO solver_used(lp);
  solver_used.SetAllAztecOptions(solver_->GetAllAztecOptions());
  solver_used.SetAllAztecParams(solver_->GetAllAztecParams());
  //
  // Set the preconditioner on the AztecOO solver object
  //
  if(Prec_.get()) {
    if(get_trace_out().get())
      trace_out() << "\nUsing a preconditioner of type \'" << typeid(*Prec_).name() << "\'\n";
    Prec_->SetUseTranspose(trans_trans(Prec_trans_,M_trans)==TRANS); // Must set mode just before use!
    solver_used.SetPrecOperator(&*Prec_);
  }
  //
  // Solve the linear system
  //
	if(get_trace_out().get())
		trace_out() << "\nSolving the linear system with AztecOO ...\n";
  int aztecStatus = solver_used.Iterate( maxIter(), relTol() );
	const int     numIterUsed = solver_used.NumIters();
	const double  tolReached  = solver_used.ScaledResidual();
	if( tolReached > relTol() ) {
		TEST_FOR_EXCEPTION(
			tolReached > minRelTol(), Solvers::Exceptions::FailureToConverge
			,"LinearOpWithSolveAztecOO::solve(...): Error, AztecOO performed numIter = "
			<< numIterUsed << " iterations (maxIter = " << maxIter() << ") and failed to converge and "
			<< "only achieved a scaled residual tolerance of " << tolReached
			<< " > minRelTol = " << minRelTol() << " > relTol = " << relTol() << "!"
			);
		if(get_trace_out().get())
			trace_out()
				<< "\nWarning! The linear solver took " << numIterUsed << " iterations "
				<< "(maxIter = " << maxIter() << ") and reached a scaled "
				<< "residual tolerance tolReached of relTol = " << relTol() << " <= tolReached = "
				<< tolReached << " <= minRelTol = " << minRelTol() << ".  Therefore, even through the "
				<< "desired tolerance was not reached we will continue anyway!\n";
	}
	else {
		if(get_trace_out().get())
			trace_out()
				<< "\nSolved the linear system(s) in " << numIterUsed << " iterations "
				<< "(maxIter = " << maxIter() << ") to a scaled "
				<< "residual tolerance of " << tolReached << " <= relTol = " << relTol() << "!\n";
	}
  //
  // Release the Epetra_MultiVector views of X and Y;
  X = Teuchos::null;
  Y = Teuchos::null;
  //
  // Scale the solution by alpha
  //
  if(alpha != 1.0 ) scale( alpha, X_inout );
  //
	if(get_trace_out().get())
		trace_out() << "\n*** Leaving LinearOpWithSolveAztecOO::solve(...) ...\n\n";
}

Teuchos::RefCountPtr<const LinearOpWithSolve<Scalar> >
LinearOpWithSolveAztecOO::clone_lows() const
{
	return Teuchos::null;  // Cloning can not fully be supported unfortunately (think about this?)
}

Teuchos::RefCountPtr<const LinearOp<Scalar> >
LinearOpWithSolveAztecOO::preconditioner() const
{
	return Teuchos::rcp(&tsfcore_Prec_,false); // This is the preconditioner!
}

} // namespace Nonlin
} // namespace TSFCore
