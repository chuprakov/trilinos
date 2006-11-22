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

#include "MoochoPack_MoochoThyraSolver.hpp"
#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_DefaultFiniteDifferenceModelEvaluator.hpp"
#include "Thyra_DefaultStateEliminationModelEvaluator.hpp"
#include "Thyra_DefaultEvaluationLoggerModelEvaluator.hpp"
#include "Thyra_SpmdMultiVectorFileIO.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace {

//
// ParameterList parameters and sublists
//

const std::string SolveMode_name = "Solve Mode";
const Teuchos::RefCountPtr<
  Teuchos::StringToIntegralParameterEntryValidator<
    MoochoPack::MoochoThyraSolver::ESolveMode
  >
>
solveModeValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<MoochoPack::MoochoThyraSolver::ESolveMode>(
    Teuchos::tuple<std::string>(
      "Forward Solve"
      ,"Optimize"
      )
    ,Teuchos::tuple<MoochoPack::MoochoThyraSolver::ESolveMode>(
      MoochoPack::MoochoThyraSolver::SOLVE_MODE_FORWARD
      ,MoochoPack::MoochoThyraSolver::SOLVE_MODE_OPTIMIZE
      )
    ,""
    )
  );
const std::string SolveMode_default = "Optimize";

const std::string NLPType_name = "NLP Type";
const Teuchos::RefCountPtr<
  Teuchos::StringToIntegralParameterEntryValidator<
    MoochoPack::MoochoThyraSolver::ENLPType
    >
  >
nlpTypeValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<MoochoPack::MoochoThyraSolver::ENLPType>(
    Teuchos::tuple<std::string>(
      "First Order"
      ,"Direct"
      )
    ,Teuchos::tuple<MoochoPack::MoochoThyraSolver::ENLPType>(
      MoochoPack::MoochoThyraSolver::NLP_TYPE_FIRST_ORDER
      ,MoochoPack::MoochoThyraSolver::NLP_TYPE_DIRECT
      )
    ,""
    )
  );
const std::string NLPType_default = "First Order";

const std::string NonlinearlyEliminateStates_name = "Nonlinearly Eliminate States";
const bool NonlinearlyEliminateStates_default = false;

const std::string UseFiniteDifferences_name = "Use Finite Differences";
const bool UseFiniteDifferences_default = false;

const std::string FiniteDifferenceSettings_name = "Finite Difference Settings";

const std::string FwdNewtonTol_name = "Forward Newton Tolerance";
const double FwdNewtonTol_default = -1.0;

const std::string FwdNewtonMaxIters_name = "Forward Newton Max Iters";
const int FwdNewtonMaxIters_default = 20;

const std::string StateGuessFileBaseName_name = "State Guess File Base Name";
const std::string StateGuessFileBaseName_default = "";

const std::string StateGuessScale_name = "State Guess Scale";
const double StateGuessScale_default = 1.0;

const std::string ParamGuessFileBaseName_name = "Parameters Guess File Base Name";
const std::string ParamGuessFileBaseName_default = "";

const std::string ParamGuessScale_name = "Parameters Guess Scale";
const double ParamGuessScale_default = 1.0;

const std::string StateSoluFileBaseName_name = "State Solution File Base Name";
const std::string StateSoluFileBaseName_default = "";

const std::string ParamSoluFileBaseName_name = "Parameters Solution File Base Name";
const std::string ParamSoluFileBaseName_default = "";

//
// Other stuff
//

typedef AbstractLinAlgPack::value_type  Scalar;

Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
readVectorFromFile(
  const std::string                                                   &fileNameBase
  ,const Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >  &vs
  ,const double                                                       scaleBy = 1.0
  )
{
  Thyra::SpmdMultiVectorFileIO<Scalar> fileIO;
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
  Thyra::SpmdMultiVectorFileIO<Scalar> fileIO;
  fileIO.writeToFile(vec,fileNameBase);
}

} // namespace

namespace MoochoPack {

// Constructors/initialization

MoochoThyraSolver::MoochoThyraSolver(
  const std::string    &paramsXmlFileName
  ,const std::string   &extraParamsXmlString
  ,const std::string   &paramsUsedXmlOutFileName
  ,const std::string   &paramsXmlFileNameOption
  ,const std::string   &extraParamsXmlStringOption
  ,const std::string   &paramsUsedXmlOutFileNameOption
  )
  :paramsXmlFileName_(paramsXmlFileName)
  ,extraParamsXmlString_(extraParamsXmlString)
  ,paramsUsedXmlOutFileName_(paramsUsedXmlOutFileName)
  ,paramsXmlFileNameOption_(paramsXmlFileNameOption)
  ,extraParamsXmlStringOption_(extraParamsXmlStringOption)
  ,paramsUsedXmlOutFileNameOption_(paramsUsedXmlOutFileNameOption)
  ,solveMode_(SOLVE_MODE_OPTIMIZE)
  ,nlpType_(NLP_TYPE_FIRST_ORDER)
  ,nonlinearlyElimiateStates_(false)
  ,use_finite_diff_(false)
  ,fwd_newton_tol_(-1.0)
  ,fwd_newton_max_iters_(20)
  ,stateGuessFileBase_("")
  ,scaleStateGuess_(1.0)
  ,paramGuessFileBase_("")
  ,scaleParamGuess_(1.0)
  ,stateSoluFileBase_("")
  ,paramSoluFileBase_("")
{}

MoochoThyraSolver::~MoochoThyraSolver()
{}

void MoochoThyraSolver::setupCLP(
  Teuchos::CommandLineProcessor *clp
  )
{
  TEST_FOR_EXCEPT(0==clp);
  solver_.setup_commandline_processor(clp);
  clp->setOption(
    paramsXmlFileNameOption().c_str(),&paramsXmlFileName_
    ,"Name of an XML file containing parameters for linear solver options to be appended first."
    );
  clp->setOption(
    extraParamsXmlStringOption().c_str(),&extraParamsXmlString_
    ,"An XML string containing linear solver parameters to be appended second."
    );
  clp->setOption(
    paramsUsedXmlOutFileNameOption().c_str(),&paramsUsedXmlOutFileName_
    ,"Name of an XML file that can be written with the parameter list after it has been used on completion of this program."
    );
}

void MoochoThyraSolver::readParameters( std::ostream *out )
{
  Teuchos::RefCountPtr<Teuchos::ParameterList>
    paramList = this->getParameterList();
  if(!paramList.get())
    paramList = Teuchos::rcp(new Teuchos::ParameterList("MoochoThyraSolver"));
  if(paramsXmlFileName().length()) {
    if(out) *out << "\nReading parameters from XML file \""<<paramsXmlFileName()<<"\" ...\n";
    Teuchos::updateParametersFromXmlFile(paramsXmlFileName(),&*paramList);
  }
  if(extraParamsXmlString().length()) {
    if(out) *out << "\nAppending extra parameters from the XML string \""<<extraParamsXmlString()<<"\" ...\n";
    Teuchos::updateParametersFromXmlString(extraParamsXmlString(),&*paramList);
  }
  if( paramsXmlFileName().length() || extraParamsXmlString().length() )
    this->setParameterList(paramList);
}

// Overridden from ParameterListAcceptor

void MoochoThyraSolver::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  if(paramList.get())
    paramList->validateParameters(*getValidParameters(),0); // Just validate my level!
  paramList_ = paramList;
  if(paramList_.get()) {
    solveMode_ = solveModeValidator->getIntegralValue(
      *paramList_,SolveMode_name,SolveMode_default);
    nlpType_ = nlpTypeValidator->getIntegralValue(
      *paramList_,NLPType_name,NLPType_default);
    nonlinearlyElimiateStates_ = paramList_->get(
      NonlinearlyEliminateStates_name,NonlinearlyEliminateStates_default);
    use_finite_diff_ = paramList_->get(
      UseFiniteDifferences_name,UseFiniteDifferences_default);
    fwd_newton_tol_ = paramList_->get(
      FwdNewtonTol_name,FwdNewtonTol_default);
    fwd_newton_max_iters_ = paramList_->get(
      FwdNewtonMaxIters_name,FwdNewtonMaxIters_default);
    stateGuessFileBase_ = paramList_->get(
      StateGuessFileBaseName_name,StateGuessFileBaseName_default);
    scaleStateGuess_ = paramList_->get(
      StateGuessScale_name,StateGuessScale_default);
    paramGuessFileBase_ = paramList_->get(
      ParamGuessFileBaseName_name,ParamGuessFileBaseName_default);
    scaleParamGuess_ = paramList_->get(
      ParamGuessScale_name,ParamGuessScale_default);
    stateSoluFileBase_ = paramList_->get(
      StateSoluFileBaseName_name,StateSoluFileBaseName_default);
    paramSoluFileBase_ = paramList_->get(
      ParamSoluFileBaseName_name,ParamSoluFileBaseName_default);
  }
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
MoochoThyraSolver::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
MoochoThyraSolver::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
MoochoThyraSolver::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
MoochoThyraSolver::getValidParameters() const
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> pl;
  if(pl.get()==NULL) {
    pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set(
      SolveMode_name,SolveMode_default
      ,"The type of solve to perform."
      ,solveModeValidator
      );
    pl->set(
      NLPType_name,NLPType_default
      ,"The type of MOOCHO NLP subclass to use."
      ,nlpTypeValidator
      );
    pl->set(
      NonlinearlyEliminateStates_name,NonlinearlyEliminateStates_default
      ,"If true, then the model's state equations and state variables\n"
      "are nonlinearlly eliminated using a forward solver."
      );
    pl->set(
      UseFiniteDifferences_name,UseFiniteDifferences_default
      ,"Use finite differences for missing model derivatives (Direct NLP only).\n"
      "See the options in the sublist \"" + FiniteDifferenceSettings_name + "\"."
      );
    if(1) {
      Teuchos::ParameterList
        &fdSublist = pl->sublist(FiniteDifferenceSettings_name);
      Thyra::DirectionalFiniteDiffCalculator<Scalar> dfdcalc;
      fdSublist.setParameters(*dfdcalc.getValidParameters());
    }
    pl->set(
      FwdNewtonTol_name,FwdNewtonTol_default
      ,"Tolarance used for the forward state solver in eliminating\n"
      "the state equations/variables."
      );
    pl->set(
      FwdNewtonMaxIters_name,FwdNewtonMaxIters_default
      ,"Maximum number of iterations allows for the forward state\n"
      "solver in eliminating the state equations/variables."
      );
    pl->set(
      StateGuessFileBaseName_name,StateGuessFileBaseName_default
      ,"The base name of a file that is used to read in the guess for\n"
      "the state variables.  If set, then a file for each process must be\n"
      "present for each process to read from."
      );
    pl->set(
      StateGuessScale_name,StateGuessScale_default
      ,"Sets the constant by which the initial guess for the state variables are scaled\n"
      "will be scaled by before they are used.  This feature allows for controlled\n"
      "experiments where a solution is perturbed and the solver must resolve the\n"
      "problem."
      );
    pl->set(
      ParamGuessFileBaseName_name,ParamGuessFileBaseName_default
      ,"The base name of a file that is used to read in the guess for\n"
      "the parameters.  If set, then a file for each process must be\n"
      "present for each process to read from."
      );
    pl->set(
      ParamGuessScale_name,ParamGuessScale_default
      ,"Sets the constant by which the initial guess for the parameters are scaled\n"
      "will be scaled by before they are used.  This feature allows for controlled\n"
      "experiments where a solution is perturbed and the solver must resolve the\n"
      "problem."
      );
    pl->set(
      StateSoluFileBaseName_name,StateSoluFileBaseName_default
      ,"If specified, a file with this basename will be written to with\n"
      "the final value of the state variables.  A different file for each\n"
      "process will be created.  Note that these files can be used for the\n"
      "initial guess for the state variables."
      );
    pl->set(
      ParamSoluFileBaseName_name,ParamSoluFileBaseName_default
      ,"If specified, a file with this basename will be written to with\n"
      "the final value of the parameters.  A different file for each\n"
      "process will be created.  Note that these files can be used for the\n"
      "initial guess for the parameters."
      );
  }
  return pl;
}

// Misc Access/Setup

void MoochoThyraSolver::setSolveMode( const ESolveMode solveMode )
{
  solveMode_ = solveMode_;
}

MoochoThyraSolver::ESolveMode
MoochoThyraSolver::getSolveMode() const
{
  return solveMode_;
}

MoochoSolver& MoochoThyraSolver::getSolver()
{
  return solver_;
}

const MoochoSolver& MoochoThyraSolver::getSolver() const
{
  return solver_;
}

// Model specification, setup, solve, and solution extraction.

void MoochoThyraSolver::setModel(
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
  //const int numProcs = Teuchos::GlobalMPISession::getNProc();

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

  switch(solveMode_) {
    case SOLVE_MODE_FORWARD: {
      nlp = rcp(new NLPFirstOrderThyraModelEvaluator(outerModel_,-1,-1));
      break;
    }
    case SOLVE_MODE_OPTIMIZE: {
      // Setup finite difference object
      RefCountPtr<Thyra::DirectionalFiniteDiffCalculator<Scalar> > direcFiniteDiffCalculator;
      if(use_finite_diff_) {
        direcFiniteDiffCalculator = rcp(new Thyra::DirectionalFiniteDiffCalculator<Scalar>());
        if(paramList_.get())
          direcFiniteDiffCalculator->setParameterList(
            Teuchos::sublist(paramList_,FiniteDifferenceSettings_name)
            );
      }
      if( nonlinearlyElimiateStates_ ) {
        // Create a Thyra::NonlinearSolverBase object to solve and eliminate the
        // state variables and the state equations
        Teuchos::RefCountPtr<Thyra::DampenedNewtonNonlinearSolver<Scalar> >
          stateSolver = rcp(new Thyra::DampenedNewtonNonlinearSolver<Scalar>()); // ToDo: Replace with MOOCHO!
        stateSolver->defaultTol(fwd_newton_tol_);
        stateSolver->defaultMaxNewtonIterations(fwd_newton_max_iters_);
        // Create the reduced Thyra::ModelEvaluator object for p -> g_hat(p)
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
        switch(nlpType_) {
          case NLP_TYPE_DIRECT: {
            Teuchos::RefCountPtr<NLPDirectThyraModelEvaluator>
              nlpDirect = rcp(new NLPDirectThyraModelEvaluator(outerModel_,p_idx,g_idx));
            if(use_finite_diff_) {
              nlpDirect->set_direcFiniteDiffCalculator(direcFiniteDiffCalculator);
            }
            nlp = nlpDirect;
            break;
          }
          case NLP_TYPE_FIRST_ORDER: {
            nlp = rcp(new NLPFirstOrderThyraModelEvaluator(outerModel_,p_idx,g_idx));
            break;
          }
          default:
            TEST_FOR_EXCEPT(true);
        }
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);
  }
    
  // Set the NLP
  solver_.set_nlp(nlp);

}

void MoochoThyraSolver::readInitialGuess(
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

void MoochoThyraSolver::setInitialGuess(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluatorBase::InArgs<value_type> > &initialGuess
  )
{
  nominalModel_->setNominalValues(initialGuess);
}

void MoochoThyraSolver::setInitialGuess(
  const Thyra::ModelEvaluatorBase::InArgs<value_type> &initialGuess
  )
{
  nominalModel_->setNominalValues(
    Teuchos::rcp(new Thyra::ModelEvaluatorBase::InArgs<value_type>(initialGuess))
    );
}
  
MoochoSolver::ESolutionStatus	MoochoThyraSolver::solve()
{
  return solver_.solve_nlp();
}

const Thyra::ModelEvaluatorBase::InArgs<value_type>&
MoochoThyraSolver::getFinalPoint() const
{
  return finalPointModel_->getFinalPoint();
}

void MoochoThyraSolver::writeFinalSolution(
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

void MoochoThyraSolver::writeParamsFile(
  const std::string &outputXmlFileName
  ) const
{
  std::string xmlOutputFile
    = ( outputXmlFileName.length() ? outputXmlFileName : paramsUsedXmlOutFileName() );
  if( paramList_.get() && xmlOutputFile.length() ) {
    Teuchos::writeParameterListToXmlFile(*paramList_,xmlOutputFile);
  }
}

} // namespace MoochoPack
