// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_RythmosSolver.hpp"
#include "Piro_ValidPiroParameters.hpp"

#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_RampingIntegrationControlStrategy.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_CompositeIntegrationObserver.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Rythmos_IntegratorBuilder.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Piro_RythmosNOX_RowSumUpdater.hpp"

#ifdef Piro_ENABLE_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

#include "Thyra_ScaledModelEvaluator.hpp"

#include <iostream>
#include <string>

template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver(
    Teuchos::RCP<Teuchos::ParameterList> appParams,
    Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > in_model,
    Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer) :
  model(in_model),
  num_p(in_model->createInArgs().Np()),
  num_g(in_model->createOutArgs().Ng()),
  out(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  *out << "\nA) Get the base parameter list ...\n";
  //

  RCP<Teuchos::ParameterList> rythmosPL = sublist(appParams, "Rythmos", true);
  rythmosPL->validateParameters(*getValidRythmosParameters(),0);

  {
    const std::string verbosity = rythmosPL->get("Verbosity Level", "VERB_DEFAULT");
    solnVerbLevel = Teuchos::VERB_DEFAULT;
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
  }

  t_final = rythmosPL->get("Final Time", 0.1);

  const std::string stepperType = rythmosPL->get("Stepper Type", "Backward Euler");

  //
  *out << "\nC) Create and initalize the forward model ...\n";
  //

  *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
  //

  if (rythmosPL->get<std::string>("Nonlinear Solver Type") == "Rythmos") {
    Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<Scalar> > rythmosTimeStepSolver =
      Rythmos::timeStepNonlinearSolver<Scalar>();
    if (rythmosPL->getEntryPtr("NonLinear Solver")) {
      RCP<Teuchos::ParameterList> nonlinePL =
	sublist(rythmosPL, "NonLinear Solver", true);
      rythmosTimeStepSolver->setParameterList(nonlinePL);
    }
    fwdTimeStepSolver = rythmosTimeStepSolver;
  }
  else if (rythmosPL->get<std::string>("Nonlinear Solver Type") == "NOX") {
#ifdef Piro_ENABLE_NOX
    Teuchos::RCP<Thyra::NOXNonlinearSolver> nox_solver =  Teuchos::rcp(new Thyra::NOXNonlinearSolver);
    Teuchos::RCP<Teuchos::ParameterList> nox_params = Teuchos::rcp(new Teuchos::ParameterList);
    *nox_params = appParams->sublist("NOX");
    nox_solver->setParameterList(nox_params);
    fwdTimeStepSolver = nox_solver;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Requested NOX solver for a Rythmos Transient solve, Trilinos was not built with NOX enabled.  Please rebuild Trilinos or use the native Rythmos nonlinear solver.");
#endif

  }

  if (stepperType == "Backward Euler") {
    fwdStateStepper = Rythmos::backwardEulerStepper<Scalar> (model, fwdTimeStepSolver);
    fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  }
  else if (stepperType == "Explicit RK") {
    fwdStateStepper = Rythmos::explicitRKStepper<Scalar>(model);
    fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  }
  else if (stepperType == "BDF") {
    Teuchos::RCP<Teuchos::ParameterList> BDFparams =
      Teuchos::sublist(rythmosPL, "Rythmos Stepper", true);
    Teuchos::RCP<Teuchos::ParameterList> BDFStepControlPL =
      Teuchos::sublist(BDFparams,"Step Control Settings");

    fwdStateStepper = Teuchos::rcp( new Rythmos::ImplicitBDFStepper<Scalar>(model,fwdTimeStepSolver,BDFparams) );
    fwdStateStepper->setInitialCondition(model->getNominalValues());

  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION( true, Teuchos::Exceptions::InvalidParameter,
				std::endl << "Error! Piro::Epetra::RythmosSolver: Invalid Steper Type: "
				<< stepperType << std::endl);

  // Step control strategy
  {
    // If the stepper can accept a step control strategy, then attempt to build one.
    RCP<Rythmos::StepControlStrategyAcceptingStepperBase<Scalar> > scsa_stepper =
      Teuchos::rcp_dynamic_cast<Rythmos::StepControlStrategyAcceptingStepperBase<Scalar> >(fwdStateStepper);

    if (Teuchos::nonnull(scsa_stepper)) {
      std::string step_control_strategy = rythmosPL->get("Step Control Strategy Type", "None");

      if (step_control_strategy == "None") {
	// don't do anything, stepper will build default
      }
      else if (step_control_strategy == "ImplicitBDFRamping") {

	const RCP<Rythmos::ImplicitBDFStepperRampingStepControl<Scalar> > rscs =
	  rcp(new Rythmos::ImplicitBDFStepperRampingStepControl<Scalar>);

	const RCP<ParameterList> p = parameterList(rythmosPL->sublist("Rythmos Step Control Strategy"));
	rscs->setParameterList(p);

	scsa_stepper->setStepControlStrategy(rscs);
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Error! Piro::Epetra::RythmosSolver: Invalid step control strategy type: " << step_control_strategy << std::endl);
      }

    }

  }

  {
    RCP<Teuchos::ParameterList>
      integrationControlPL = sublist(rythmosPL, "Rythmos Integration Control", true);

    RCP<Rythmos::DefaultIntegrator<Scalar> > defaultIntegrator;

    if (rythmosPL->get("Rythmos Integration Control Strategy", "Simple") == "Simple") {
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::simpleIntegrationControlStrategy<Scalar>(integrationControlPL));
    }
    else if(rythmosPL->get<std::string>("Rythmos Integration Control Strategy") == "Ramping") {
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::rampingIntegrationControlStrategy<Scalar>(integrationControlPL));
    }

    fwdStateIntegrator = defaultIntegrator;
  }

  fwdStateIntegrator->setParameterList(sublist(rythmosPL, "Rythmos Integrator", true));

  if (Teuchos::nonnull(observer)) {
    fwdStateIntegrator->setIntegrationObserver(observer);
  }
}


template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver(
    const Teuchos::RCP<Rythmos::DefaultIntegrator<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Rythmos::StepperBase<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<Scalar> > &underlyingModel,
    Scalar finalTime,
    const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<Scalar> > &icModel,
    Teuchos::EVerbosityLevel verbosityLevel) :
  fwdStateIntegrator(stateIntegrator),
  fwdStateStepper(stateStepper),
  fwdTimeStepSolver(timeStepSolver),
  model(underlyingModel),
  initialConditionModel(icModel),
  num_p(model->createInArgs().Np()),
  num_g(model->createOutArgs().Ng()),
  t_final(finalTime),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel(verbosityLevel)
{
  fwdStateStepper->setNonconstModel(underlyingModel);
}


template <typename Scalar>
Teuchos::RCP<const Rythmos::IntegratorBase<Scalar> >
Piro::RythmosSolver<Scalar>::getRythmosIntegrator() const
{
  return fwdStateIntegrator;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::RythmosSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::RythmosSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::RythmosSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::RythmosSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_space(j);
  else return model->get_x_space(); // j == num_g
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::RythmosSolver<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> Piro::RythmosSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> Piro::RythmosSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g + 1);

  const int l = 0;
  // Solution sensitivity
  outArgs.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
      num_g,
      l,
      Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

  return outArgs;
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // TODO: Support more than 1 parameter and 1 response
  const int j = 0;
  const int l = 0;

  // Parse InArgs
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(l);
  }

  // Parse OutArgs
  RCP<Thyra::VectorBase<Scalar> > g_out;
  if (num_g > 0) {
    g_out = outArgs.get_g(j);
  }
  const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);

  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = model->getNominalValues();

  // Set paramters p_in as part of initial conditions
  if (num_p > 0){
    state_ic.set_p(l, p_in);
  }

  *out << "\nstate_ic:\n" << Teuchos::describe(state_ic,solnVerbLevel);

  RCP<Thyra::MultiVectorBase<Scalar> > dgxdp_out;
  const Thyra::ModelEvaluatorBase::DerivativeSupport dgxdp_support =
    outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, num_g, l);
  if (dgxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
    const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgxdp_deriv =
      outArgs.get_DgDp(num_g, l);
    dgxdp_out = dgxdp_deriv.getMultiVector();
  }

  const bool requestedSensitivities = Teuchos::nonnull(dgxdp_out);

  RCP<const Thyra::VectorBase<Scalar> > finalSolution;
  if (!requestedSensitivities) {
    //
    *out << "\nE) Solve the forward problem ...\n";
    //

    fwdStateStepper->setInitialCondition(state_ic);
    fwdStateIntegrator->setStepper(fwdStateStepper, t_final, true);

    Teuchos::Array<RCP<const Thyra::VectorBase<Scalar> > > x_final_array;
    fwdStateIntegrator->getFwdPoints(
        Teuchos::tuple<Scalar>(t_final), &x_final_array, NULL, NULL);
    finalSolution = x_final_array[0];

    if (Teuchos::VERB_MEDIUM <= solnVerbLevel) {
      std::cout << "Final Solution\n" << *finalSolution << std::endl;
    }

    // As post-processing step, calculate responses at final solution
    Thyra::ModelEvaluatorBase::InArgs<Scalar> model_inargs = model->createInArgs();
    model_inargs.set_x(finalSolution);
    if (num_p > 0) {
      model_inargs.set_p(l, p_in);
    }

    Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
    if (Teuchos::nonnull(g_out)) {
      Thyra::put_scalar(Teuchos::ScalarTraits<Scalar>::zero(), g_out.ptr());
      model_outargs.set_g(j, g_out);
    }

    model->evalModel(model_inargs, model_outargs);

  } else { // Computing sensitivities
    //
    *out << "\nE) Solve the forward problem with Sensitivities...\n";
    //
    RCP<Rythmos::ForwardSensitivityStepper<Scalar> > stateAndSensStepper =
      Rythmos::forwardSensitivityStepper<Scalar>();
    stateAndSensStepper->initializeSyncedSteppers(
        model, l, model->getNominalValues(),
        fwdStateStepper, fwdTimeStepSolver);

    //
    // Set the initial condition for the state and forward sensitivities
    //

    const RCP< Thyra::VectorBase<Scalar> > s_bar_init =
      createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
    const RCP< Thyra::VectorBase<Scalar> > s_bar_dot_init =
      createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());

    if (Teuchos::is_null(initialConditionModel)) {
      // The initial condition is assumed to be independent from the parameters
      // Therefore, the initial condition for the sensitivity is zero
      assign(s_bar_init.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
    } else {
      // Use initialConditionModel to compute initial condition for sensitivity
      Thyra::ModelEvaluatorBase::InArgs<Scalar> initCondInArgs = initialConditionModel->createInArgs();
      initCondInArgs.set_p(l, inArgs.get_p(l));

      Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs = initialConditionModel->createOutArgs();
      typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
      const RCP<DMVPV> s_bar_init_downcasted = Teuchos::rcp_dynamic_cast<DMVPV>(s_bar_init);
      const Thyra::ModelEvaluatorBase::Derivative<Scalar> initCond_deriv(
          s_bar_init_downcasted->getNonconstMultiVector(),
          Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
      initCondOutArgs.set_DgDp(initCondOutArgs.Ng() - 1, l, initCond_deriv);

      initialConditionModel->evalModel(initCondInArgs, initCondOutArgs);
    }
    assign(s_bar_dot_init.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

    RCP<const Rythmos::StateAndForwardSensitivityModelEvaluator<Scalar> >
      stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

    Thyra::ModelEvaluatorBase::InArgs<Scalar>
      state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

    // Copy time, parameters etc.
    state_and_sens_ic.setArgs(state_ic);
    // Set initial condition for x_bar = [ x; s_bar ]
    state_and_sens_ic.set_x(
        stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
        );
    // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
    state_and_sens_ic.set_x_dot(
        stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(), s_bar_dot_init)
        );

    stateAndSensStepper->setInitialCondition(state_and_sens_ic);

    //
    // Use a StepperAsModelEvaluator to integrate the state+sens
    //

    const RCP<Rythmos::StepperAsModelEvaluator<Scalar> >
      stateAndSensIntegratorAsModel = Rythmos::stepperAsModelEvaluator(
          Teuchos::rcp_implicit_cast<Rythmos::StepperBase<Scalar> >(stateAndSensStepper),
          Teuchos::rcp_implicit_cast<Rythmos::IntegratorBase<Scalar> >(fwdStateIntegrator),
          state_and_sens_ic);
    const int solutionResponseIndex = 0;

    *out << "\nUse the StepperAsModelEvaluator to integrate state + sens x_bar(p,t_final) ... \n";
    Teuchos::OSTab tab(out);

    // Solution sensitivity in column-oriented (Jacobian) MultiVector form
    RCP<const Thyra::MultiVectorBase<Scalar> > dxdp;

    const RCP<Thyra::VectorBase<Scalar> > x_bar_final =
      Thyra::createMember(stateAndSensIntegratorAsModel->get_g_space(j));
    // Extract pieces of x_bar_final to prepare output
    {
      const RCP<const Thyra::ProductVectorBase<Scalar> > x_bar_final_downcasted =
        Thyra::productVectorBase<Scalar>(x_bar_final);

      // Solution
      const int solutionBlockIndex = 0;
      finalSolution = x_bar_final_downcasted->getVectorBlock(solutionBlockIndex);

      // Sensitivity
      const int sensitivityBlockIndex = 1;
      const RCP<const Thyra::VectorBase<Scalar> > s_bar_final =
        x_bar_final_downcasted->getVectorBlock(sensitivityBlockIndex);
      {
        typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
        const RCP<const DMVPV> s_bar_final_downcasted = Teuchos::rcp_dynamic_cast<const DMVPV>(s_bar_final);

        dxdp = s_bar_final_downcasted->getMultiVector();
      }
    }

    eval_g(
        *stateAndSensIntegratorAsModel,
        l, *state_ic.get_p(l),
        t_final,
        solutionResponseIndex, &*x_bar_final
        );

    *out
      << "\nx_bar_final = x_bar(p,t_final) evaluated using "
      << "stateAndSensIntegratorAsModel:\n"
      << Teuchos::describe(*x_bar_final,solnVerbLevel);

    if (Teuchos::nonnull(dgxdp_out)) {
      Thyra::assign(dgxdp_out.ptr(), *dxdp);
    }

    *out << "\nF) Check the solution to the forward problem ...\n";
  }

  // Return the final solution as an additional g-vector, if requested
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*finalSolution, gx_out.ptr());
  }
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::RythmosSolver<Scalar>::getValidRythmosParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidRythmosParams"));

  validPL->set<std::string>("Nonlinear Solver Type", "");

  Teuchos::setStringToIntegralParameter<int>(
    "Nonlinear Solver Type",
    "Rythmos",
    "Determines which nonlinear solver to use.",
    Teuchos::tuple<std::string>("Rythmos","NOX"),
    validPL.get()
    );

  validPL->sublist("NonLinear Solver", false, "");
  validPL->sublist("Rythmos Builder", false, "");


  validPL->set<double>("Final Time", 1.0, "");
  validPL->sublist("Rythmos Stepper", false, "");
  validPL->sublist("Rythmos Integrator", false, "");
  validPL->set<std::string>("Rythmos Integration Control Strategy", "Simple", "");
  validPL->set<std::string>("Step Control Strategy Type", "None", "");
  validPL->sublist("Rythmos Step Control Strategy", false, "");
  validPL->sublist("Rythmos Integration Control", false, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<std::string>("Stepper Type", "", "");

  validPL->set<double>("Alpha", 1.0, "");
  validPL->set<double>("Beta", 1.0, "");
  validPL->set<double>("Max State Error", 1.0, "");
  validPL->set<std::string>("Name", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");

  return validPL;
}
