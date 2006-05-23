// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "MoochoPack_ThyraModelEvaluatorSolver.hpp"
#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_DefaultFiniteDifferenceModelEvaluator.hpp"
#include "Thyra_DefaultStateEliminationModelEvaluator.hpp"
#include "Thyra_DefaultEvaluationLoggerModelEvaluator.hpp"
#include "Thyra_ParallelMultiVectorFileIO.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace {

typedef AbstractLinAlgPack::value_type  Scalar;

Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
readVectorFromFile(
  const std::string                                                   &fileNameBase
  ,const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >  &vs
  ,const double                                                       scaleBy = 1.0
  )
{
  Thyra::ParallelMultiVectorFileIO<Scalar> fileIO;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    vec = fileIO.readVectorFromFile(fileNameBase,vs);
  Thyra::Vt_S(&*vec,scaleBy);
  return vec;
}

void writeVectorToFile(
  const Thyra::VectorBase<Scalar>    &vec
  ,const std::string                 &fileNameBase
  )
{
  Thyra::ParallelMultiVectorFileIO<Scalar> fileIO;
  fileIO.writeToFile(vec,fileNameBase);
}

const int numFDMethodTypes = 7;

Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType fdMethodTypeValues[numFDMethodTypes] =
{
  Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_ONE
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_CENTRAL
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_AUTO
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_CENTRAL
  ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO
};

const char* fdMethodTypeNames[numFDMethodTypes] =
{
  "order-one"
  ,"order-two"
  ,"order-two-central"
  ,"order-two-auto"
  ,"order-four"
  ,"order-four-central"
  ,"order-four-auto"
};

} // namespace


namespace MoochoPack {

ThyraModelEvaluatorSolver::ThyraModelEvaluatorSolver()
  :do_sim_(false)
  ,use_direct_(false)
  ,use_black_box_(false)
  ,use_finite_diff_(false)
  ,fd_method_type_(Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_ONE)
  ,fd_step_len_(1e-4)
  ,fwd_newton_tol_(-1.0)
  ,fwd_newton_max_iters_(20)
  ,stateGuessFileBase_("")
  ,scaleStateGuess_(1.0)
  ,paramGuessFileBase_("")
  ,scaleParamGuess_(1.0)
  ,stateSoluFileBase_("")
  ,paramSoluFileBase_("")
{}

void ThyraModelEvaluatorSolver::setupCLP(
  Teuchos::CommandLineProcessor *clp
  )
{
  solver_.setup_commandline_processor(clp);
  clp->setOption( "do-sim", "do-opt",  &do_sim_, "Flag for if only the square constraints are solved" );
  clp->setOption( "use-direct", "use-first-order",  &use_direct_, "Flag for if we use the NLPDirect or NLPFirstOrderInfo implementation." );
  clp->setOption( "black-box", "all-at-once",  &use_black_box_, "Flag for if we do black-box NAND or all-at-once SAND." );
  clp->setOption( "use-finite-diff", "no-use-finite-diff",  &use_finite_diff_, "Flag for if we are using finite differences or not." );
  clp->setOption(
    "finite-diff-type", &fd_method_type_
    ,numFDMethodTypes,fdMethodTypeValues,fdMethodTypeNames
    ,"The type of finite differences used [---use-finite-diff]"
    );
  clp->setOption( "fd-step-len", &fd_step_len_, "The finite-difference step length ( < 0 determine automatically )." );
  clp->setOption( "fwd-newton-tol", &fwd_newton_tol_, "Nonlinear residual tolerance for the forward Newton solve." );
  clp->setOption( "fwd-newton-max-iters", &fwd_newton_max_iters_, "Max number of iterations for the forward Newton solve." );
  clp->setOption( "x-guess-file", &stateGuessFileBase_, "Base file name to read the guess of the state x from." );
  clp->setOption( "scale-x-guess", &scaleStateGuess_, "Amount to scale the guess for x read in by --x-guess-file." );
  clp->setOption( "p-guess-file", &paramGuessFileBase_, "Base file name to read the guess of the parameters p from." );
  clp->setOption( "scale-p-guess", &scaleParamGuess_, "Amount to scale the guess for p read in by --p-guess-file." );
  clp->setOption( "x-solu-file", &stateSoluFileBase_, "Base file name to write the state solution x to." );
  clp->setOption( "p-solu-file", &paramSoluFileBase_, "Base file name to write the parameter solution p to." );
}

MoochoSolver& ThyraModelEvaluatorSolver::getSolver()
{
  return solver_;
}

const MoochoSolver& ThyraModelEvaluatorSolver::getSolver() const
{
  return solver_;
}

void ThyraModelEvaluatorSolver::setModel(
  const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> > &model
  ,const int                                                     p_idx
  ,const int                                                     g_idx
  )
{

  using Teuchos::rcp;
  using Teuchos::RefCountPtr;
  using NLPInterfacePack::NLP;
  using NLPInterfacePack::NLPDirectThyraModelEvaluator;
  using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;

  origModel_ = model;
  p_idx_ = p_idx;
  g_idx_ = g_idx;

  const int procRank = Teuchos::GlobalMPISession::getRank();
  const int numProcs = Teuchos::GlobalMPISession::getNProc();

  //
  // Wrap the orginal model
  //

  Teuchos::RefCountPtr<std::ostream>
    modelEvalLogOut;
  if(procRank==0)
    modelEvalLogOut = rcp(new std::ofstream("ModelEvaluationLog.out"));
  else
    modelEvalLogOut = rcp(new Teuchos::oblackholestream());
  Teuchos::RefCountPtr<Thyra::DefaultEvaluationLoggerModelEvaluator<Scalar> >
    loggerThyraModel
    = rcp(
      new Thyra::DefaultEvaluationLoggerModelEvaluator<Scalar>(
        origModel_,modelEvalLogOut
        )
      );
  nominalModel_
    = rcp(
      new Thyra::DefaultNominalBoundsOverrideModelEvaluator<Scalar>(loggerThyraModel,Teuchos::null)
      );
  finalPointModel_
    = rcp(
      new Thyra::DefaultFinalPointCaptureModelEvaluator<value_type>(nominalModel_)
      );
  outerModel_ = finalPointModel_;

  //
  // Create the NLP
  //
    
  Teuchos::RefCountPtr<NLP> nlp;

  if(do_sim_) {
    nlp = rcp(new NLPFirstOrderThyraModelEvaluator(outerModel_,-1,-1));
  }
  else {
    // Setup finite difference object
    RefCountPtr<Thyra::DirectionalFiniteDiffCalculator<Scalar> > direcFiniteDiffCalculator;
    if(use_finite_diff_) {
      direcFiniteDiffCalculator = rcp(new Thyra::DirectionalFiniteDiffCalculator<Scalar>());
      direcFiniteDiffCalculator->fd_method_type(fd_method_type_);
      direcFiniteDiffCalculator->fd_step_size(fd_step_len_);
    }
    if( use_black_box_ ) {
      // Create a Thyra::NonlinearSolverBase object to solve and eliminate the
      // state variables and the state equations
      Teuchos::RefCountPtr<Thyra::DampenedNewtonNonlinearSolver<Scalar> >
        stateSolver = rcp(new Thyra::DampenedNewtonNonlinearSolver<Scalar>()); // ToDo: Replace with MOOCHO!
      stateSolver->defaultTol(fwd_newton_tol_);
      stateSolver->defaultMaxNewtonIterations(fwd_newton_max_iters_);
      // Create the reduced  Thyra::ModelEvaluator object for p -> g_hat(p)
      Teuchos::RefCountPtr<Thyra::DefaultStateEliminationModelEvaluator<Scalar> >
        reducedThyraModel = rcp(new Thyra::DefaultStateEliminationModelEvaluator<Scalar>(outerModel_,stateSolver));
      Teuchos::RefCountPtr<Thyra::ModelEvaluator<Scalar> >
        finalReducedThyraModel;
      if(use_finite_diff_) {
        // Create the finite-difference wrapped Thyra::ModelEvaluator object
        Teuchos::RefCountPtr<Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar> >
          fdReducedThyraModel = rcp(
            new Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar>(
              reducedThyraModel,direcFiniteDiffCalculator)
            );
        finalReducedThyraModel = fdReducedThyraModel;
      }
      else {
        finalReducedThyraModel = reducedThyraModel;
      }
      // Wrap the reduced NAND Thyra::ModelEvaluator object in an NLP object
      nlp = rcp(new NLPFirstOrderThyraModelEvaluator(finalReducedThyraModel,p_idx,g_idx));
    }
    else {
      if(use_direct_) {
        Teuchos::RefCountPtr<NLPDirectThyraModelEvaluator>
          nlpDirect = rcp(new NLPDirectThyraModelEvaluator(outerModel_,p_idx,g_idx));
        if(use_finite_diff_) {
          nlpDirect->set_direcFiniteDiffCalculator(direcFiniteDiffCalculator);
        }
        nlp = nlpDirect;
      }
      else {
        nlp = rcp(new NLPFirstOrderThyraModelEvaluator(outerModel_,p_idx,g_idx));
      }
    }
  }
    
  // Set the NLP
  solver_.set_nlp(nlp);

}

void ThyraModelEvaluatorSolver::readInitialGuess(
  std::ostream *out_arg
  )
{
  using Teuchos::OSTab;
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  Teuchos::RefCountPtr<Thyra::ModelEvaluatorBase::InArgs<value_type> >
    initialGuess = Teuchos::rcp(new Thyra::ModelEvaluatorBase::InArgs<value_type>(origModel_->createInArgs()));
  if(stateGuessFileBase_ != "") {
    if(out.get())
      *out << "\nReading the guess of the state \'x\' from the file(s) with base name \""<<stateGuessFileBase_<<"\" ...\n";
    initialGuess->set_x(readVectorFromFile(stateGuessFileBase_,origModel_->get_x_space(),scaleStateGuess_));
  }
  else {
    if(out.get())
      *out << "\nNot reading the guess of the state \'x\' from a file!\n";
  }
  if(paramGuessFileBase_ != "") {
    if(out.get())
      *out << "\nReading the guess of the parameters \'p\' from the file(s) with base name \""<<paramGuessFileBase_<<"\" ...\n";
    initialGuess->set_p(p_idx_,readVectorFromFile(paramGuessFileBase_,origModel_->get_p_space(p_idx_),scaleParamGuess_));
  }
  else {
    if(out.get())
      *out << "\nNot reading the guess of the parameters \'p\' from a file!\n";
  }
  nominalModel_->setNominalValues(initialGuess);
}
  
MoochoSolver::ESolutionStatus	ThyraModelEvaluatorSolver::solve()
{
  return solver_.solve_nlp();
}

void ThyraModelEvaluatorSolver::writeFinalSolution(
  std::ostream *out_arg
  ) const
{
  using Teuchos::OSTab;
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  if( stateSoluFileBase_ != "" && finalPointModel_->getFinalPoint().get_x().get() ) {
    if(out.get())
      *out << "\nWriting the state solution \'x\' to the file(s) with base name \""<<stateSoluFileBase_<<"\" ...\n";
    writeVectorToFile(*finalPointModel_->getFinalPoint().get_x(),stateSoluFileBase_);
  }
  if( paramSoluFileBase_ != "" ) {
    if(out.get())
      *out << "\nWriting the parameter solution \'p\' to the file(s) with base name \""<<paramSoluFileBase_<<"\" ...\n";
    writeVectorToFile(*finalPointModel_->getFinalPoint().get_p(p_idx_),paramSoluFileBase_);
  }
}

} // namespace MoochoPack
