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
// TSFCoreNonlinLinearOpWithSolveAztecOO.cpp

#include "TSFCoreNonlinLinearOpWithSolveAztecOO.hpp"
#include "TSFCore_get_Epetra_MultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Time.hpp"

namespace {

inline
Teuchos::ETransp convert( TSFCore::ETransp trans_in )
{
	Teuchos::ETransp  trans_out;
	switch(trans_in) {
		case TSFCore::NOTRANS:
			trans_out = Teuchos::NO_TRANS;
			break;
		case TSFCore::TRANS:
			trans_out = Teuchos::TRANS;
			break;
		default:
			assert(0); // Should never get here!
	}
	return trans_out;
}

} // namespace

namespace TSFCore {
namespace Nonlin {

// Constructors / initializers / accessors

LinearOpWithSolveAztecOO::LinearOpWithSolveAztecOO(
	const int                    maxIter
	,const double                relTol
	,const double                minRelTol
	)
	:maxIter_(maxIter)
	,relTol_(relTol)
	,minRelTol_(minRelTol)
{
	resetCounters();
}

LinearOpWithSolveAztecOO::~LinearOpWithSolveAztecOO()
{
	if(get_trace_out().get())
		trace_out()
			<< "\nLinearOpWithSolveAztecOO::~LinearOpWithSolveAztecOO(): Printing final stats:"
			<< "\n  Number of forward solves           = " << numFwdSolves_
			<< "\n  Number of forward solve iterations = " << numFwdLinearIters_
			<< "\n  forward solve seconds              = " << fwdLinearCPU_
			<< "\n  Number of adjoint solves           = " << numAdjSolves_
			<< "\n  Number of adjoint solve iterations = " << numAdjLinearIters_
			<< "\n  adjoint solve seconds              = " << adjLinearCPU_
			<< std::endl;
}

void LinearOpWithSolveAztecOO::initialize(
	const Teuchos::RefCountPtr<Epetra_Operator>      &Op
	,const ETransp                                   Op_trans
	,const Teuchos::RefCountPtr<AztecOO>             &solver
	,const Teuchos::RefCountPtr<Epetra_Operator>     &Prec
	,const ETransp                                   Prec_trans
	,const Epetra::ProductOperator::EApplyMode       Prec_inverse
	,const bool                                      adjointSupported
	)
{
	Teuchos::Time timer("");
#ifdef _DEBUG
	const char func_name[] = "LinearOpWithSolveAztecOO::initialize(...)";
	TEST_FOR_EXCEPTION(Op.get()==NULL,std::invalid_argument,func_name<<": Error!");
	TEST_FOR_EXCEPTION(solver.get()==NULL,std::invalid_argument,func_name<<": Error!");
  // ToDo: Check vector spaces!
#endif
  if(get_trace_out().get()) {
    trace_out()
			<< "\nLinearOpWithSolveAztecOO::initialize(...):"
			<< "\n  Using an operator of type \'" << typeid(*Op).name() << "\'";
		if(Prec.get())
			trace_out()
				<< "\nLinearOpWithSolveAztecOO::initialize(...):"
				<< "\n  Using an preconditioner of type \'" << typeid(*Prec).name() << "\'";
    trace_out() << std::endl;
	}
	// Set the references
  Op_ = Op;
  Op_trans_ = Op_trans;
  tsfcore_Op_.initialize(Op,Op_trans);
  solver_ = solver;
  if(Prec.get()) tsfcore_Prec_.initialize(Prec,Prec_trans);
  Prec_ = Prec;
  Prec_trans_ = Prec_trans;
	Prec_inverse_ = Prec_inverse;
  adjointSupported_ = adjointSupported;
	// Scale the operator and/or the preconditioner
	Teuchos::RefCountPtr<const Epetra::LinearSystemScaler::Scaling> scaling;
/*
	if( Epetra_RowMatrix *Op_rm = dynamic_cast<Epetra_RowMatrix*>(&*Op_) ) {
		if(get_trace_out().get())
			trace_out()
				<< "\nLinearOpWithSolveAztecOO::initialize(...): "
				<< "Operator support Epetra_RowMatrix, computing row and column scaling ...";
		timer.start(true);
		scaling	= linearSystemScaler().computeScaling(*Op_rm,convert(Op_trans_));
		timer.stop();
		if(get_trace_out().get()) {
			double left_norm[1], right_norm[1];
			scaling->leftScaling()->NormInf(left_norm);
			scaling->rightScaling()->NormInf(right_norm);
			trace_out()
				<< "\n  => time = " << timer.totalElapsedTime() << " sec\n"
				<< "\n  ||leftScaling||inf  = " << left_norm[0]
				<< "\n  ||rightScaling||inf = " << right_norm[0]
				<< std::endl;
		}
	}
	else {
*/
		scaling = Teuchos::rcp(new Epetra::LinearSystemScaler::Scaling()); // No scaling!
/*
	}
*/
	// Create scaled objects for the operator and the preconditioner.
	// Note, AztecOO calls ApplyInverse(...)  to apply the
	// preconditioner when using an external preconditioner.
	if(get_trace_out().get())
		trace_out()
			<< "\nLinearOpWithSolveAztecOO::initialize(...): "
			<< "Generating the scaled operator and preconditioner for the forward solve  ...\n";
	linearSystemScaler().generateSolveOps(
		*scaling, convert(Op_trans_)
		,Op_, convert(Op_trans_)
		,Prec_, convert(Prec_trans_), Prec_inverse_
		,Teuchos::NO_TRANS
		,&fwd_Op_
		,Teuchos::RIGHT_SIDE,Epetra::ProductOperator::APPLY_MODE_APPLY_INVERSE, &fwd_Prec_
		,get_trace_out().get()
		);
	if(adjointSupported_) {
		if(get_trace_out().get())
			trace_out()
				<< "\nLinearOpWithSolveAztecOO::initialize(...): "
				<< "Generating the scaled operator and preconditioner for the adjoint solve  ...\n";
		linearSystemScaler().generateSolveOps(
			*scaling, convert(Op_trans_)
			,Op_, convert(Op_trans_)
			,Prec_, convert(Prec_trans_), Prec_inverse_
			,Teuchos::TRANS
			,&adj_Op_
			,Teuchos::RIGHT_SIDE,Epetra::ProductOperator::APPLY_MODE_APPLY_INVERSE, &adj_Prec_
			,get_trace_out().get()
			);
	}
}
	
void LinearOpWithSolveAztecOO::setUninitialized(
	Teuchos::RefCountPtr<Epetra_Operator>            *Op
	,ETransp                                         *Op_trans
	,Teuchos::RefCountPtr<AztecOO>                   *solver
	,Teuchos::RefCountPtr<Epetra_Operator>           *Prec
	,ETransp                                         *Prec_trans
	,Epetra::ProductOperator::EApplyMode             *Prec_inverse
	,bool                                            *adjointSupported
  )
{
  Op_ = Teuchos::null;
  Op_trans_ = NOTRANS;
  tsfcore_Op_.setUninitialized(Op,Op_trans);
  if(solver) *solver = solver_; solver_ = Teuchos::null;
  Prec_ = Teuchos::null;
  Prec_trans_ = NOTRANS;
  Prec_inverse_ = Epetra::ProductOperator::APPLY_MODE_APPLY;
  tsfcore_Prec_.setUninitialized(Prec,Prec_trans);
  if(adjointSupported) *adjointSupported = adjointSupported_;
}

void LinearOpWithSolveAztecOO::resetCounters()
{
	numFwdSolves_ = numFwdLinearIters_ = numAdjSolves_ = numAdjLinearIters_ = 0;
	fwdLinearCPU_ = adjLinearCPU_ = 0.0;
}

// Overridden from OpBase

Teuchos::RefCountPtr<const VectorSpace<LinearOpWithSolveAztecOO::Scalar> >
LinearOpWithSolveAztecOO::domain() const
{
  return tsfcore_Op_.domain();
}

Teuchos::RefCountPtr<const VectorSpace<LinearOpWithSolveAztecOO::Scalar> >
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
#ifdef TSFCORE_VECTOR_DERIVE_FROM_MULTI_VECTOR
	solve(M_trans,static_cast<const MultiVector<Scalar>&>(y),static_cast<MultiVector<Scalar>*>(x),1.0,convTester);
#else
	const MultiVectorCols<Scalar>  Y(Teuchos::rcp(const_cast<Vector<Scalar>*>(&y),false));
	MultiVectorCols<Scalar>        X(Teuchos::rcp(x,false));
	solve(M_trans,Y,&X,1.0,convTester);
#endif
}

void LinearOpWithSolveAztecOO::solve(
	const ETransp                         M_trans
	,const MultiVector<Scalar>            &Y_in
	,MultiVector<Scalar>                  *X_inout
	,const Scalar                         alpha
	,Solvers::ConvergenceTester<Scalar>   *convTester
	) const
{
	Teuchos::Time timer("");
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
	// Get the operator and preconditioner being used for this solve
	//
	Teuchos::RefCountPtr<Epetra_Operator>
		using_Op   = ( M_trans==NOTRANS ? fwd_Op_   : adj_Op_   ),
		using_Prec = ( M_trans==NOTRANS ? fwd_Prec_ : adj_Prec_ );
	//
	// Scale the RHS
	//
	Epetra_MultiVector Y_scaled = *Y; // ToDo: Make data member?
	linearSystemScaler().preSolveTransformRhs(*using_Op,convert(M_trans),&Y_scaled);
  //
  // Create the linear problem
  //
  if(get_trace_out().get())
    trace_out() << "\nUsing an operator of type \'" << typeid(*using_Op).name() << "\'\n";
  //Op_->SetUseTranspose(trans_trans(Op_trans_,M_trans)==TRANS); // Must set mode just before use!
  using_Op->SetUseTranspose(false); // Must set mode just before use!
  Epetra_LinearProblem
    lp(
      &*using_Op
      ,&*X
      ,&Y_scaled
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
  if( using_Prec.get() ) {
    if(get_trace_out().get())
      trace_out() << "\nUsing a preconditioner of type \'" << typeid(*using_Prec).name() << "\'\n";
    //Prec_->SetUseTranspose(trans_trans(Prec_trans_,M_trans)==TRANS); // Must set mode just before use!
    using_Prec->SetUseTranspose(false); // Must set mode just before use!
    solver_used.SetPrecOperator( &*using_Prec );
  }
	else {
    if(get_trace_out().get())
      trace_out() << "\nNo preconditioner specified!\n";
	}
  //
  // Solve the linear system
  //
	if(get_trace_out().get())
		trace_out()
			<< "\nSolving the linear system (M_trans = "
			<< (M_trans==NOTRANS ? "NOTRANS" : "TRANS") << ") with AztecOO ...\n";
	timer.start(true);
	solver_used.Iterate( maxIter(), relTol() ); // We ignore the returned status (see below)
	timer.stop();
	if(get_trace_out().get())
		trace_out()
			<< "\n  => time = " << timer.totalElapsedTime() << " sec\n";
	const int     numIterUsed = solver_used.NumIters();
	const double  tolReached  = solver_used.ScaledResidual();
	if(M_trans==NOTRANS) {
		++numFwdSolves_;
		numFwdLinearIters_ += numIterUsed;
		fwdLinearCPU_ += timer.totalElapsedTime();
	}
	else {
		++numAdjSolves_;
		numAdjLinearIters_ += numIterUsed;
		adjLinearCPU_ += timer.totalElapsedTime();
	}
	const double *AZ_status = solver_used.GetAztecStatus();
	if(get_trace_out().get()) {
		if(AZ_status[AZ_why]==AZ_normal) {
			trace_out() << "\nAztec returned AZ_normal\n";
		}
		else if(AZ_status[AZ_why]==AZ_param) {
			trace_out() << "\nAztec returned AZ_param\n";
		}
		else if(AZ_status[AZ_why]==AZ_breakdown) {
			trace_out() << "\nAztec returned AZ_breakdown\n";
		}
		else if(AZ_status[AZ_why]==AZ_loss) {
			trace_out() << "\nAztec returned AZ_loss\n";
		}
		else if(AZ_status[AZ_why]==AZ_ill_cond) {
			trace_out() << "\nAztec returned AZ_ill_cond\n";
		}
		else if(AZ_status[AZ_why]==AZ_maxits) {
			trace_out() << "\nAztec returned AZ_maxits\n";
		}
		else {
			trace_out() << "\nAztec returned an unknown status?\n";
		}
	}
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
	// ToDo: Act on bad Aztec return status values
	//
	// Scale the solution
	//
	linearSystemScaler().postSolveTransformSolu(
		*using_Op
		,convert(M_trans)
		,&*X
		);
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

Teuchos::RefCountPtr<const LinearOpWithSolve<LinearOpWithSolveAztecOO::Scalar> >
LinearOpWithSolveAztecOO::clone_lows() const
{
	return Teuchos::null;  // Cloning can not fully be supported unfortunately (think about this?)
}

Teuchos::RefCountPtr<const LinearOp<LinearOpWithSolveAztecOO::Scalar> >
LinearOpWithSolveAztecOO::preconditioner() const
{
	return Teuchos::rcp(&tsfcore_Prec_,false); // This is the preconditioner!
}

} // namespace Nonlin
} // namespace TSFCore
