// /////////////////////////////////////////////////////////////////////////
// rSQPAlgo_ConfigMamaJama.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <assert.h>

#include <sstream>
#include <typeinfo>
#include <iostream>

#include "debug.hpp"

#include "rSQPAlgo_ConfigMamaJama.hpp"
#include "ReducedSpaceSQPPack/src/rSQPAlgo.hpp"
#include "ReducedSpaceSQPPack/src/rSQPAlgoContainer.hpp"
#include "SparseLinAlgPack/src/MatrixSymPosDefCholFactor.hpp"                     // rHL 
//#include "ConstrainedOptimizationPack/src/MatrixSymPosDefInvCholFactor.hpp"		// .
#include "ConstrainedOptimizationPack/src/MatrixSymPosDefLBFGS.hpp"				// .
//#include "ConstrainedOptimizationPack/src/MatrixHessianSuperBasicInitDiagonal.hpp"// | rHL (super basics)
//#include "SparseLinAlgPack/src/MatrixSymDiagStd.hpp"                          // |

#include "ConstrainedOptimizationPack/src/VariableBoundsTester.hpp"

#include "NLPInterfacePack/src/NLPFirstOrderDirect.hpp"

#include "NLPInterfacePack/src/CalcFiniteDiffProd.hpp"
#include "NLPInterfacePack/src/NLPVarReductPerm.hpp"

// line search
#include "ConstrainedOptimizationPack/src/DirectLineSearchArmQuad_Strategy.hpp"
#include "ConstrainedOptimizationPack/src/DirectLineSearchArmQuad_StrategySetOptions.hpp"
#include "ConstrainedOptimizationPack/src/MeritFuncNLPL1.hpp"
#include "ConstrainedOptimizationPack/src/MeritFuncNLPModL1.hpp"

// Basis permutations and direct sparse solvers
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
#include "ConstrainedOptimizationPack/src/DecompositionSystemVarReductPerm.hpp"
#endif

#include "ConstrainedOptimizationPack/src/QPSolverRelaxedTester.hpp"
#include "ConstrainedOptimizationPack/src/QPSolverRelaxedTesterSetOptions.hpp"
#include "ConstrainedOptimizationPack/src/QPSolverRelaxedQPSchur.hpp"
#include "ConstrainedOptimizationPack/src/QPSolverRelaxedQPSchurSetOptions.hpp"
#include "ConstrainedOptimizationPack/src/QPSchurInitKKTSystemHessianFull.hpp"
//#include "ConstrainedOptimizationPack/src/QPSchurInitKKTSystemHessianSuperBasic.hpp"
#include "ConstrainedOptimizationPack/src/QPSolverRelaxedQPKWIK.hpp"
#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
//#include "ConstrainedOptimizationPack/src/QPSolverRelaxedQPOPT.hpp"
#endif

#include "ReducedSpaceSQPPack/src/std/rSQPAlgorithmStepNames.hpp"

/*
#include "ReducedSpaceSQPPack/src/std/EvalNewPointStd_StepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/EvalNewPointTailoredApproach_StepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/EvalNewPointTailoredApproachCoordinate_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/EvalNewPointTailoredApproachOrthogonal_Step.hpp"
*/

#include "ReducedSpaceSQPPack/src/std/ReducedGradientStd_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/InitFinDiffReducedHessian_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/InitFinDiffReducedHessian_StepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/ReducedHessianSecantUpdateStd_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/ReducedHessianSecantUpdateBFGSFull_Strategy.hpp"


//#include "ReducedSpaceSQPPack/src/std/ReducedHessianSecantUpdateBFGSProjected_Strategy.hpp"
//#include "ReducedSpaceSQPPack/src/std/ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.hpp"
//#include "ReducedSpaceSQPPack/src/std/ReducedHessianSecantUpdateLPBFGS_Strategy.hpp"
//#include "ReducedSpaceSQPPack/src/std/ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/BFGSUpdate_Strategy.hpp"
#include "ReducedSpaceSQPPack/src/std/BFGSUpdate_StrategySetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/RangeSpaceStepStd_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/CheckDescentRangeSpaceStep_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/CheckDecompositionFromPy_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/CheckDecompositionFromRPy_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/NullSpaceStepWithoutBounds_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/NullSpaceStepWithInequStd_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/NullSpaceStepWithInequStd_StepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/SetDBoundsStd_AddedStep.hpp"
#include "ReducedSpaceSQPPack/src/std/QPFailureReinitReducedHessian_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/CalcDFromYPYZPZ_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/CalcDFromZPZ_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/LineSearchFailureNewDecompositionSelection_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/LineSearchFilter_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/LineSearchFilter_StepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/LineSearchFullStep_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/LineSearchDirect_Step.hpp"
//#include "ReducedSpaceSQPPack/src/std/LineSearch2ndOrderCorrect_Step.hpp"
//#include "ReducedSpaceSQPPack/src/std/LineSearch2ndOrderCorrect_StepSetOptions.hpp"
//#include "ReducedSpaceSQPPack/src/std/FeasibilityStepReducedStd_Strategy.hpp"
//#include "ReducedSpaceSQPPack/src/std/FeasibilityStepReducedStd_StrategySetOptions.hpp"
//#include "ReducedSpaceSQPPack/src/std/QuasiRangeSpaceStepStd_Strategy.hpp"
//#include "ReducedSpaceSQPPack/src/std/QuasiRangeSpaceStepTailoredApproach_Strategy.hpp"
//#include "ReducedSpaceSQPPack/src/std/LineSearchWatchDog_Step.hpp"
//#include "ReducedSpaceSQPPack/src/std/LineSearchWatchDog_StepSetOptions.hpp"
//#include "ReducedSpaceSQPPack/src/std/LineSearchFullStepAfterKIter_Step.hpp"
//#include "ReducedSpaceSQPPack/src/std/CalcLambdaIndepStd_AddedStep.hpp"
#include "ReducedSpaceSQPPack/src/std/CalcReducedGradLagrangianStd_AddedStep.hpp"
#include "ReducedSpaceSQPPack/src/std/CheckConvergenceStd_AddedStep.hpp"
#include "ReducedSpaceSQPPack/src/std/CheckConvergenceStd_Strategy.hpp"
#include "ReducedSpaceSQPPack/src/std/CheckSkipBFGSUpdateStd_StepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/MeritFunc_DummyUpdate_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.hpp"
#include "ReducedSpaceSQPPack/src/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.hpp"
//#include "ReducedSpaceSQPPack/src/std/MeritFunc_PenaltyParamUpdateWithMult_AddedStep.hpp"
//#include "ReducedSpaceSQPPack/src/std/MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.hpp"
//#include "ReducedSpaceSQPPack/src/std/MeritFunc_ModifiedL1LargerSteps_AddedStep.hpp"
//#include "ReducedSpaceSQPPack/src/std/MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.hpp"
//#include "ReducedSpaceSQPPack/src/std/ActSetStats_AddedStep.hpp"
//#include "ReducedSpaceSQPPack/src/std/NumFixedDepIndep_AddedStep.hpp"

#include "ReducedSpaceSQPPack/src/std/act_set_stats.hpp"
#include "ReducedSpaceSQPPack/src/std/qp_solver_stats.hpp"
#include "ReducedSpaceSQPPack/src/std/quasi_newton_stats.hpp"

//#include "SparseLinAlgPack/src/sparse_bounds.hpp"

// Misc utilities
#include "AbstractFactoryStd.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "ThrowException.hpp"

// Stuff to read in options
#include "StringToIntMap.hpp"
#include "StringToBool.hpp"

// Stuff for exact reduced hessian
//#include "ReducedSpaceSQPPack/src/std/ReducedHessianExactStd_Step.hpp"
//#include "ReducedSpaceSQPPack/src/std/CrossTermExactStd_Step.hpp"
//#include "ReducedSpaceSQPPack/src/std/DampenCrossTermStd_Step.hpp"

namespace {
	const double INF_BASIS_COND_CHANGE_FRAC      = 1e+20;
}

namespace ReducedSpaceSQPPack {

//
// Here is where we define the default values for the algorithm.  These
// should agree with what are in the rSQPpp.opt.rSQPAlgo_ConfigMamaJama file.
//
rSQPAlgo_ConfigMamaJama::SOptionValues::SOptionValues()
	:max_basis_cond_change_frac_(-1.0)
	,exact_reduced_hessian_(false)
	,quasi_newton_(QN_AUTO)
	,num_lbfgs_updates_stored_(-1)
	,lbfgs_auto_scaling_(true)
	,hessian_initialization_(INIT_HESS_AUTO)
	,qp_solver_type_(QP_AUTO)
	,reinit_hessian_on_qp_fail_(true)
	,line_search_method_(LINE_SEARCH_AUTO)
	,merit_function_type_(MERIT_FUNC_AUTO)
	,l1_penalty_param_update_(L1_PENALTY_PARAM_AUTO)
	,full_steps_after_k_(-1)
{}

rSQPAlgo_ConfigMamaJama::rSQPAlgo_ConfigMamaJama()
{}

rSQPAlgo_ConfigMamaJama::~rSQPAlgo_ConfigMamaJama()
{}

// overridden from rSQPAlgo_Config

void rSQPAlgo_ConfigMamaJama::set_options( const options_ptr_t& options )
{
	options_ = options;
	decomp_sys_step_builder_.set_options(options);
}

const rSQPAlgo_Config::options_ptr_t&
rSQPAlgo_ConfigMamaJama::get_options() const
{
	return options_;
}

void rSQPAlgo_ConfigMamaJama::config_algo_cntr(
	rSQPAlgoContainer   *algo_cntr
	,std::ostream       *trase_out
	)
{
	namespace afp = MemMngPack;
	namespace mmp = MemMngPack;
	using mmp::ref_count_ptr;
	using DynamicCastHelperPack::dyn_cast;

	if(trase_out) {
		*trase_out
			<< std::endl
			<< "*****************************************************************\n"
			<< "*** rSQPAlgo_ConfigMamaJama configuration                     ***\n"
			<< "***                                                           ***\n"
			<< "*** Here, summary information about how the algorithm is      ***\n"
			<< "*** configured is printed so that the user can see how the    ***\n"
			<< "*** properties of the NLP and the set options influence       ***\n"
			<< "*** how an algorithm is configured.                           ***\n"
			<< "*****************************************************************\n";
	}

	// ////////////////////////////////////////////////////////////
	// A. ???

	// /////////////////////////////////////////////////////////////////////////
	// B. Create an algo object, give to algo_cntr, then give algo_cntr to algo

	if(trase_out)
		*trase_out << "\n*** Creating the rSQPAlgo algo object ...\n";

	typedef mmp::ref_count_ptr<rSQPAlgo>	algo_ptr_t;
	algo_ptr_t algo = mmp::rcp(new rSQPAlgo);
	assert(algo.get());
	algo_cntr->set_algo(algo);
	algo->set_algo_cntr(algo_cntr);

	// /////////////////////////////////////////////
	// C. Configure algo

	// /////////////////////////////////////////////////////
	// C.0 Set the nlp and track objects

	if(trase_out)
		*trase_out << "\n*** Setting the NLP and track objects to the algo object ...\n";

	algo->set_nlp( algo_cntr->get_nlp().get() );
	algo->set_track( algo_cntr->get_track() );

	// ////////////////////////////////////////////////
	// Determine what the options are:

	// Readin the options
	if(options_.get()) {
		readin_options( *options_, &uov_, trase_out );
	}
	else {
		if(trase_out) {
			*trase_out
				<< "\n*** Warning, no OptionsFromStream object was set so a default set"
					" of options will be used!\n";
		}
	}	

	NLP &nlp = algo->nlp();
	nlp.initialize(algo->algo_cntr().check_results());
	// Get the dimensions of the NLP
	const size_type
		n   = nlp.n(),
		m   = nlp.m(),
		mI  = nlp.mI(),
		r   = m, // ToDo: Compute this for real!
		dof = n - r,
		nb  = nlp.num_bounded_x();

	// Process the NLP
	NLPFirstOrderInfo    *nlp_foi = NULL;
	NLPSecondOrderInfo   *nlp_soi = NULL;
	NLPFirstOrderDirect  *nlp_fod = NULL;
	bool                 tailored_approach = false;
	decomp_sys_step_builder_.process_nlp_and_options(
		trase_out, nlp
		,&nlp_foi, &nlp_soi, &nlp_fod, &tailored_approach
		);

	const int max_dof_quasi_newton_dense
		= decomp_sys_step_builder_.current_option_values().max_dof_quasi_newton_dense_;

	// Make sure that we can handle this type of NLP currently
	THROW_EXCEPTION(
		n == m, std::logic_error
		,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
		"can not currently solve a square system of equations!" );
	THROW_EXCEPTION(
		mI, std::logic_error
		,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
		"can not currently solve an NLP with general inequalities!" );
	
	// //////////////////////////////////////////////////////
	// C.1.  Sort out the options

	if(trase_out)
		*trase_out
			<< "\n*** Sorting out some of the options given input options ...\n";

	if( m ) {
		if( tailored_approach ) {
			// Change the options for the tailored approach. 
			if(trase_out) {
				*trase_out
					<< "\nThis is a tailored approach NLP (NLPFirstOrderDirect) which forces the following options:\n"
					<< "merit_function_type         = L1;\n"
					<< "l1_penalty_parameter_update = MULT_FREE;\n"
					<< "null_space_matrix           = EXPLICIT;\n"
					;
			}
			cov_.merit_function_type_
				= MERIT_FUNC_L1;
			cov_.l1_penalty_param_update_
				= L1_PENALTY_PARAM_MULT_FREE;
			decomp_sys_step_builder_.current_option_values().null_space_matrix_type_
				= DecompositionSystemStateStepBuilderStd::NULL_SPACE_MATRIX_EXPLICIT;
		}
		
		if( !tailored_approach && uov_.merit_function_type_ != MERIT_FUNC_L1  ) {
			if(trase_out) {
				*trase_out
					<< "\nThe only merit function currently supported is:\n"
					<< "merit_function_type         = L1;\n"
					;
			}
			cov_.merit_function_type_		= MERIT_FUNC_L1;
		}
	}
	else {
		if(trase_out) {
			*trase_out
				<< "\nThere are now equality constraints (m == 0) which forces the following options:\n"
				<< "line_search_method          = DIRECT;\n"
				<< "merit_function_type         = L1;\n"
				;
		}
		cov_.line_search_method_       = LINE_SEARCH_DIRECT;
		cov_.merit_function_type_      = MERIT_FUNC_L1;
	}

	// Decide what type of quasi-newton update to use
	switch( uov_.quasi_newton_ ) {
		case QN_AUTO: {
			if(trase_out)
				*trase_out
					<< "\nquasi_newton == AUTO:"
					<< "\nnlp.num_bounded_x() == " << nlp.num_bounded_x() << ":\n";
			if( n - r > max_dof_quasi_newton_dense ) {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " > max_dof_quasi_newton_dense = "
						<< max_dof_quasi_newton_dense <<  ":\n"
						<< "setting quasi_newton == LBFGS\n";
				cov_.quasi_newton_ = QN_LBFGS;
			}
			else {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " <= max_dof_quasi_newton_dense = "
						<< max_dof_quasi_newton_dense << ":\n"
						<< "setting quasi_newton == BFGS\n";
				cov_.quasi_newton_ = QN_BFGS;
			}
			break;
		}
		case QN_BFGS:
		case QN_PBFGS:
		case QN_LBFGS:
		case QN_LPBFGS:
			cov_.quasi_newton_ = uov_.quasi_newton_;
			break;
	    default:
			assert(0); // Invalid option!
	}

	if( uov_.qp_solver_type_ == QP_AUTO && nb == 0 && mI == 0 ) {
		cov_.qp_solver_type_ = QP_AUTO;
	}

	// ToDo: Sort out the rest of the options!

	// Set the default options that where not already set yet
	set_default_options(uov_,&cov_,trase_out);

	// ToDo: Implement the 2nd order correction linesearch
	if( cov_.line_search_method_ == LINE_SEARCH_2ND_ORDER_CORRECT ) {
		if(trase_out)
			*trase_out <<
				"\nline_search_method == 2ND_ORDER_CORRECT:\n"
				"Sorry, the second order corrrection linesearch is not updated yet!\n"
				"setting line_search_method = DIRECT ...\n";
		cov_.line_search_method_ = LINE_SEARCH_DIRECT;
	}
	if( cov_.line_search_method_ == LINE_SEARCH_WATCHDOG ) {
		if(trase_out)
			*trase_out <<
				"\nline_search_method ==WATCHDOG:\n"
				"Sorry, the watchdog linesearch is not updated yet!\n"
				"setting line_search_method = DIRECT ...\n";
		cov_.line_search_method_ = LINE_SEARCH_DIRECT;
	}
	
	// Adjust the Quasi-Newton if QPKWIK is used!
	if( cov_.qp_solver_type_ == QP_QPKWIK && ( nb || mI ) ) {
		if(trase_out)
			*trase_out
				<< "\nqp_solver == QPKWIK and nlp.num_bounded_x() == " << nb << " and nlp.mI() == " << mI << ":\n"
				<< "Setting quasi_newton == BFGS...\n";
		cov_.quasi_newton_ = QN_BFGS;
	}

	// /////////////////////////////////////////////////////
	// C.1. Create the decomposition system object

	typedef ref_count_ptr<DecompositionSystem> decomp_sys_ptr_t;
	decomp_sys_ptr_t decomp_sys;
	decomp_sys_step_builder_.create_decomp_sys(
		trase_out, nlp, nlp_foi, nlp_soi, nlp_fod, tailored_approach
		,&decomp_sys
		);

#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
	ref_count_ptr<DecompositionSystemVarReductPerm>
		decomp_sys_perm = mmp::rcp_dynamic_cast<DecompositionSystemVarReductPerm>(decomp_sys);
#endif

	// /////////////////////////////////////////////////////
	// C.2. Create and set the state object

	if(trase_out)
		*trase_out
			<< "\n*** Creating the state object and setting up iteration quantity objects ...\n";

	{
		//
		// Create the state object with the vector spaces
		//

		typedef ref_count_ptr<rSQPState>   state_ptr_t;
		state_ptr_t
			state = mmp::rcp(
				new rSQPState(
					decomp_sys
					,nlp.space_x()
					,nlp.space_c()
					,nlp.space_h()
					,( m
					   ? ( tailored_approach
						   ? ( nlp_fod->var_dep().size() 
							   ? nlp.space_x()->sub_space(nlp_fod->var_dep())->clone()
							   : mmp::null )
						   : decomp_sys->space_range() // could be NULL for BasisSystemPerm
						   )
					   : mmp::null
						)
					,( m
					   ? ( tailored_approach
						   ?( nlp_fod->var_indep().size()
							  ? nlp.space_x()->sub_space(nlp_fod->var_indep())->clone()
							  : mmp::null )
						   : decomp_sys->space_null() // could be NULL for BasisSystemPerm
						   )
					   : nlp.space_x()
						)
					)
				);
				
		//
		// Set the iteration quantities for the NLP matrix objects
		//

		decomp_sys_step_builder_.add_iter_quantities(
			trase_out, nlp, nlp_foi, nlp_soi, nlp_fod, tailored_approach, decomp_sys
			,state
			);

		// Add reduced Hessian

		if( !cov_.exact_reduced_hessian_ ) {
			ref_count_ptr<afp::AbstractFactory<MatrixSymOp> >
				abstract_factory_rHL = mmp::null;
			// Only maintain the orginal matrix if we have inequality constraints and therefore will be
			// needing a QP solver (which may be QPSchur which needs an accurate original matrix for
			// iterative refinement).
			const bool
				maintain_original = ( nb || mI );
			// Maintain the factor if a QP solver is needed, QPSchur is being used, or we are checking
			// results
			const bool
				maintain_inverse = ( (!nb && !mI && m==r) || cov_.qp_solver_type_==QP_QPSCHUR
									 || algo->algo_cntr().check_results() );
			switch( cov_.quasi_newton_ ) {
				case QN_BFGS:
					abstract_factory_rHL = mmp::rcp(
						new afp::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor,MatrixSymPosDefCholFactor::PostMod>(
							MatrixSymPosDefCholFactor::PostMod(
								maintain_original      // maintain_original
								,maintain_inverse      // maintain_factor
								,true                  // allow_factor      (always!)
								)
							)
						);
					break;
				case QN_LBFGS:
					abstract_factory_rHL = mmp::rcp(
						new afp::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefLBFGS,MatrixSymPosDefLBFGS::PostMod>(
							MatrixSymPosDefLBFGS::PostMod(
								cov_.num_lbfgs_updates_stored_  //
								,maintain_original              // maintain_original
								,maintain_inverse               // maintain_inverse
								,cov_.lbfgs_auto_scaling_       // 
								)
							)
						);
					break;
				default:
					assert(0); // Should not be called for now!
			}
			
			state->set_iter_quant(
				rHL_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<MatrixSymOp>(
						1
						,rHL_name
						,abstract_factory_rHL
						)
					)
				);
		}
		else {
			assert(0); // ToDo: Add rHL for an exact reduced Hessian!
		}
		
		//
		// Set the NLP merit function 
		//

		if( cov_.line_search_method_ != LINE_SEARCH_NONE 
			 && cov_.line_search_method_ != LINE_SEARCH_FILTER) {
			ref_count_ptr<afp::AbstractFactory<MeritFuncNLP> >
				merit_func_factory = mmp::null;
			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1:
					merit_func_factory = mmp::rcp(
						new afp::AbstractFactoryStd<MeritFuncNLP,MeritFuncNLPL1>());
					break;
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
					merit_func_factory = mmp::rcp(
						new afp::AbstractFactoryStd<MeritFuncNLP,MeritFuncNLPModL1>());
					break;
				default:
					assert(0);	// local programming error
			}
			state->set_iter_quant(
				merit_func_nlp_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<MeritFuncNLP>(
						1
						,merit_func_nlp_name
						,merit_func_factory
						)
					)
				);
		}

		if (cov_.line_search_method_ == LINE_SEARCH_FILTER)
		  {
		    // Add the filter iteration quantity
		    state->set_iter_quant(
				FILTER_IQ_STRING
				,mmp::rcp(
				     new IterQuantityAccessContiguous<Filter_T>(1,FILTER_IQ_STRING)
				)
		    );
		  }

		if( nb || mI ) {
			// Add bounds on d
			state->set_iter_quant(
				dl_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<VectorMutable>(
						2, dl_name
						, nlp.space_x() ) )
				);
			state->set_iter_quant(
				du_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<VectorMutable>(
						2, du_name
						, nlp.space_x() ) )
				);
			// Add active-set iteration quantity
			state->set_iter_quant(
				act_set_stats_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<ActSetStats>( 1, act_set_stats_name ) )
				);
			// Add QP solver stats iteration quantity
			state->set_iter_quant(
				qp_solver_stats_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<QPSolverStats>( 1, qp_solver_stats_name ) )
				);
		}

		//
		// Resize the number of storage locations (these can be changed later).
		//
		// Also, touch all of the value_type, index_type and vector iteration quantities
		// that we know about so that when state.dump_iter_quant() is called, all of the
		// iteration quantities will be included.
		//
		
		typedef IterQuantityAccessContiguous<value_type>            IQ_scalar_cngs;
		typedef IterQuantityAccessContiguous<VectorMutable>   IQ_vector_cngs;

		dyn_cast<IQ_vector_cngs>(state->x()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->f()).resize(2);
		if(m) dyn_cast<IQ_vector_cngs>(state->c()).resize(2);
		state->Gf();
		if(mI) dyn_cast<IQ_vector_cngs>(state->h()).resize(2);
		if(m && nlp_foi) state->Gc();
		if(mI) state->Gh();

		if( m
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
		   && decomp_sys_perm.get() == NULL
#endif
			) state->py();
		if(m) dyn_cast<IQ_vector_cngs>(state->Ypy()).resize(2);
		if( m
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
			&& decomp_sys_perm.get() == NULL
#endif
			) state->pz();
		if(m) dyn_cast<IQ_vector_cngs>(state->Zpz()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->d()).resize(2);

		if( n > m ) {
			dyn_cast<IQ_vector_cngs>(state->rGf()).resize(2);
			state->w();
			state->zeta();
			state->qp_grad();
		}
		state->eta();
		
		dyn_cast<IQ_scalar_cngs>(state->alpha()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->mu()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->phi()).resize(2);

		dyn_cast<IQ_scalar_cngs>(state->opt_kkt_err()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->feas_kkt_err()).resize(2);
		if( n > m ) {
			dyn_cast<IQ_vector_cngs>(state->rGL()).resize(2);
		}
		if(m) dyn_cast<IQ_vector_cngs>(state->lambda()).resize(2);
		if(mI) dyn_cast<IQ_vector_cngs>(state->lambdaI()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->nu()).resize(2);

		// Set the state object
		algo->set_state( state );
	}

	// /////////////////////////////////////////////////////
	// C.3  Create and set the step objects

	if(trase_out)
		*trase_out << "\n*** Creating and setting the step objects ...\n";

//	typedef ref_count_ptr<QPSolverRelaxed> qp_solver_ptr_t;
//	qp_solver_ptr_t qp_solver;
//	typedef ConstrainedOptimizationPack::QPSolverRelaxedTester QPSolverRelaxedTester;
//	typedef ref_count_ptr<QPSolverRelaxedTester> qp_tester_ptr_t;
//	qp_tester_ptr_t qp_tester;
//	typedef ref_count_ptr<FeasibilityStepReducedStd_Strategy> feasibility_step_strategy_ptr_t;
//	feasibility_step_strategy_ptr_t  feasibility_step_strategy;

	{

		//
		// Create some standard step objects that will be used by many different
		// specific algorithms
		//
		
		typedef ref_count_ptr<AlgorithmStep>   algo_step_ptr_t;

		// Create the EvalNewPoint step and associated objects
		algo_step_ptr_t                                    eval_new_point_step           = mmp::null;
		ref_count_ptr<CalcFiniteDiffProd>                  calc_fd_prod                  = mmp::null;
		ref_count_ptr<VariableBoundsTester>                bounds_tester                 = mmp::null;
		ref_count_ptr<NewDecompositionSelection_Strategy>  new_decomp_selection_strategy = mmp::null;
		decomp_sys_step_builder_.create_eval_new_point(
			trase_out, nlp, nlp_foi, nlp_soi, nlp_fod, tailored_approach, decomp_sys
			,&eval_new_point_step, &calc_fd_prod, &bounds_tester, &new_decomp_selection_strategy
			);

		// RangeSpace_Step
		algo_step_ptr_t    range_space_step_step = mmp::null;
		if( !tailored_approach ) {
			range_space_step_step = mmp::rcp(new RangeSpaceStepStd_Step());
		}

		// Check and change decomposition 
		algo_step_ptr_t    check_decomp_from_py_step  = mmp::null;
		algo_step_ptr_t    check_decomp_from_Rpy_step = mmp::null;
		if( new_decomp_selection_strategy.get() && cov_.max_basis_cond_change_frac_ < INF_BASIS_COND_CHANGE_FRAC ) {
			check_decomp_from_py_step = mmp::rcp(
				new CheckDecompositionFromPy_Step(
					new_decomp_selection_strategy
					,cov_.max_basis_cond_change_frac_
					) );
			check_decomp_from_Rpy_step = mmp::rcp(
				new CheckDecompositionFromRPy_Step(
					new_decomp_selection_strategy
					,cov_.max_basis_cond_change_frac_
					) );
		}

		// CheckDescentRangeSpaceStep
		algo_step_ptr_t    check_descent_range_space_step_step = mmp::null;
		if( algo->algo_cntr().check_results() ) {
			check_descent_range_space_step_step = mmp::rcp(new CheckDescentRangeSpaceStep_Step(calc_fd_prod));
		}

		// ReducedGradient_Step
		algo_step_ptr_t    reduced_gradient_step = mmp::null;
		if( !tailored_approach ) {
			reduced_gradient_step = mmp::rcp(new ReducedGradientStd_Step());
		}

		// CheckSkipBFGSUpdate
		algo_step_ptr_t    check_skip_bfgs_update_step = mmp::null;
		if(!cov_.exact_reduced_hessian_) {
			ref_count_ptr<CheckSkipBFGSUpdateStd_Step>
				step = mmp::rcp(new CheckSkipBFGSUpdateStd_Step());
			if(options_.get()) {
				CheckSkipBFGSUpdateStd_StepSetOptions
					opt_setter( step.get() );
				opt_setter.set_options( *options_ );
			}
			check_skip_bfgs_update_step = step;
		}

		// ReducedHessian_Step
		algo_step_ptr_t    reduced_hessian_step = mmp::null;
		{
			// Get the strategy object that will perform the actual secant update.
			ref_count_ptr<ReducedHessianSecantUpdate_Strategy>
				secant_update_strategy = mmp::null;
			switch( cov_.quasi_newton_ )
			{
			    case QN_BFGS:
			    case QN_PBFGS:
			    case QN_LBFGS:
			    case QN_LPBFGS:
				{
					// create and setup the actual BFGS strategy object
					typedef ref_count_ptr<BFGSUpdate_Strategy> bfgs_strategy_ptr_t;
					bfgs_strategy_ptr_t
						bfgs_strategy = mmp::rcp(new BFGSUpdate_Strategy);
					if(options_.get()) { 
						BFGSUpdate_StrategySetOptions
							opt_setter( bfgs_strategy.get() );
						opt_setter.set_options( *options_ );
					}
					switch( cov_.quasi_newton_ ) {
					    case QN_BFGS:
					    case QN_LBFGS:
						{
							secant_update_strategy = mmp::rcp(new ReducedHessianSecantUpdateBFGSFull_Strategy(bfgs_strategy));
							break;
						}
					    case QN_PBFGS:
					    case QN_LPBFGS:
						{
							THROW_EXCEPTION(
								true, std::logic_error
								,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
								"The quansi_newton options of PBFGS and LPBFGS have not been updated yet!" );
							break;
						}
					}
					break;
				}
			    default:
					assert(0);
			}
			
			// Finally build the step object
			reduced_hessian_step = mmp::rcp(
				new ReducedHessianSecantUpdateStd_Step( secant_update_strategy ) );
			// Add the QuasiNewtonStats iteration quantity
			algo->state().set_iter_quant(
				quasi_newton_stats_name
				,mmp::rcp(new IterQuantityAccessContiguous<QuasiNewtonStats>(
					1
					,quasi_newton_stats_name
#ifdef _MIPS_CXX
					,mmp::ref_count_ptr<MemMngPack::AbstractFactoryStd<QuasiNewtonStats,QuasiNewtonStats> >(
						new MemMngPack::AbstractFactoryStd<QuasiNewtonStats,QuasiNewtonStats>())
#endif
					)
				));
		}

		// InitFinDiffReducedHessian_Step
		algo_step_ptr_t  init_red_hess_step = mmp::null;
		if( cov_.hessian_initialization_ != INIT_HESS_IDENTITY ) {
				Algorithm::poss_type poss;
			assert(poss = algo->get_step_poss( ReducedHessian_name ) );
			InitFinDiffReducedHessian_Step::EInitializationMethod
				init_hess;
			switch( cov_.hessian_initialization_ ) {
				case INIT_HESS_FIN_DIFF_SCALE_IDENTITY:
					init_hess = InitFinDiffReducedHessian_Step::SCALE_IDENTITY;
					break;
				case INIT_HESS_FIN_DIFF_SCALE_DIAGONAL:
					init_hess = InitFinDiffReducedHessian_Step::SCALE_DIAGONAL;
					break;
				case INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS:
					init_hess = InitFinDiffReducedHessian_Step::SCALE_DIAGONAL_ABS;
					break;
				default:
					assert(0);	// only local programming error?
			}
			mmp::ref_count_ptr<InitFinDiffReducedHessian_Step>
				_init_red_hess_step = mmp::rcp(new InitFinDiffReducedHessian_Step(init_hess));
			if(options_.get()) {
				InitFinDiffReducedHessian_StepSetOptions
					opt_setter( _init_red_hess_step.get() );
				opt_setter.set_options( *options_ );
			}
			init_red_hess_step = _init_red_hess_step;
		}

		// NullSpace_Step
		algo_step_ptr_t    set_d_bounds_step    = mmp::null;
		algo_step_ptr_t    null_space_step_step = mmp::null;
		if( mI == 0 && nb == 0 ) {
			null_space_step_step = mmp::rcp(new NullSpaceStepWithoutBounds_Step());
		}
		else {
			// Step object that sets bounds for QP subproblem
			set_d_bounds_step = mmp::rcp(new SetDBoundsStd_AddedStep());
			// QP Solver object
			mmp::ref_count_ptr<QPSolverRelaxed>  qp_solver = mmp::null;
			switch( cov_.qp_solver_type_ ) {
				case QP_QPSCHUR: {
					using ConstrainedOptimizationPack::QPSolverRelaxedQPSchur;
					using ConstrainedOptimizationPack::QPSolverRelaxedQPSchurSetOptions;
					using ConstrainedOptimizationPack::QPSchurInitKKTSystemHessianFull;
//					using ConstrainedOptimizationPack::QPSchurInitKKTSystemHessianSuperBasic;
					mmp::ref_count_ptr<ConstrainedOptimizationPack::QPSolverRelaxedQPSchur::InitKKTSystem>
						init_kkt_sys = mmp::null;
					switch( cov_.quasi_newton_ ) {
						case QN_BFGS:
						case QN_LBFGS:
							init_kkt_sys = mmp::rcp(new QPSchurInitKKTSystemHessianFull());
							break;
						case QN_PBFGS:
						case QN_LPBFGS:
							THROW_EXCEPTION(true,std::logic_error,"Error! PBFGS and LPBFGS are not updated yet!");
/*
							init_kkt_sys = new QPSchurInitKKTSystemHessianSuperBasic();
*/
							break;
						default:
							assert(0);
					}
					mmp::ref_count_ptr<QPSolverRelaxedQPSchur>
						_qp_solver = mmp::rcp(new QPSolverRelaxedQPSchur(init_kkt_sys));
					QPSolverRelaxedQPSchurSetOptions
						qp_options_setter(_qp_solver.get());
					qp_options_setter.set_options( *options_ );
					qp_solver = _qp_solver; // give ownership to delete!
					break;
				}
				case QP_QPKWIK: {
					using ConstrainedOptimizationPack::QPSolverRelaxedQPKWIK;
					mmp::ref_count_ptr<QPSolverRelaxedQPKWIK>
						_qp_solver = mmp::rcp(new QPSolverRelaxedQPKWIK());
					qp_solver = _qp_solver; // give ownership to delete!
					break;
				}
				case QP_QPOPT: {
					THROW_EXCEPTION(true,std::logic_error,"Error! QPKWIK interface is not updated yet!");
/*
#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
					using ConstrainedOptimizationPack::QPSolverRelaxedQPOPT;
					QPSolverRelaxedQPOPT
						*_qp_solver = new QPSolverRelaxedQPOPT();
					qp_solver = _qp_solver; // give ownership to delete!
#else
					throw std::logic_error(
						"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : QPOPT is not supported,"
						" must define CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT!" );
#endif
*/
					break;
				}
				default:
					assert(0);
			}
			// QP solver tester
			mmp::ref_count_ptr<QPSolverRelaxedTester> 
				qp_solver_tester = mmp::rcp(new QPSolverRelaxedTester());
			if(options_.get()) {
				QPSolverRelaxedTesterSetOptions
					opt_setter( qp_solver_tester.get() );
				opt_setter.set_options( *options_ );
			}
			// The null-space step
			mmp::ref_count_ptr<NullSpaceStepWithInequStd_Step>
				null_space_step_with_inequ_step = mmp::rcp(
					new NullSpaceStepWithInequStd_Step(
						qp_solver, qp_solver_tester ) );
			if(options_.get()) {
				NullSpaceStepWithInequStd_StepSetOptions
					opt_setter( null_space_step_with_inequ_step.get() );
				opt_setter.set_options( *options_ );
			}
			null_space_step_step = null_space_step_with_inequ_step;
			// Step for reinitialization reduced Hessian on QP failure
			null_space_step_step = mmp::rcp(
				new QPFailureReinitReducedHessian_Step(null_space_step_step)
				);
		}

		// CalcDFromYPYZPZ_Step
		algo_step_ptr_t    calc_d_from_Ypy_Zpy_step = mmp::null;
		{
			calc_d_from_Ypy_Zpy_step = mmp::rcp(new CalcDFromYPYZPZ_Step());
		}

		// CalcReducedGradLagrangianStd_AddedStep
		algo_step_ptr_t    calc_reduced_grad_lagr_step = mmp::null;
		{
			calc_reduced_grad_lagr_step = mmp::rcp(
				new CalcReducedGradLagrangianStd_AddedStep() );
		}

		// CheckConvergence_Step
		algo_step_ptr_t    check_convergence_step = mmp::null;
			{
			// Create the strategy object
			ref_count_ptr<CheckConvergenceStd_Strategy>
				check_convergence_strategy = mmp::rcp(new CheckConvergenceStd_Strategy());

			if(options_.get()) 
				{
				CheckConvergence_StrategySetOptions
					opt_setter( check_convergence_strategy.get() );
				opt_setter.set_options( *options_ );
				}
			
			ref_count_ptr<CheckConvergenceStd_AddedStep>
				_check_convergence_step = mmp::rcp(new CheckConvergenceStd_AddedStep(check_convergence_strategy));
			
			check_convergence_step = _check_convergence_step;
			}

		// MeritFuncPenaltyParamUpdate_Step
		algo_step_ptr_t    merit_func_penalty_param_update_step = mmp::null;
		if( cov_.line_search_method_ == LINE_SEARCH_FILTER ) {
		  // We don't need to update a penalty parameter for the filter method :-)
		}
		else if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<MeritFunc_PenaltyParamUpdate_AddedStep>
				param_update_step = mmp::null;
			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1: {
					switch(cov_.l1_penalty_param_update_) {
						case L1_PENALTY_PARAM_WITH_MULT:
//							param_update_step
//								= mmp::rcp(new  MeritFunc_PenaltyParamUpdateWithMult_AddedStep());
							THROW_EXCEPTION(
								true, std::logic_error
								,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
								"The l1_penalty_parameter_update option of MULT_FREE has not been updated yet!" );
							break;
						case L1_PENALTY_PARAM_MULT_FREE:
							param_update_step
								= mmp::rcp(new  MeritFunc_PenaltyParamUpdateMultFree_AddedStep());
							break;
						default:
							assert(0);
					}
					break;
				}
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
//					param_update_step = new  MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(
//											mmp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
					THROW_EXCEPTION(
						true, std::logic_error
						,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
						"The merit_function_type options of MODIFIED_L1 and MODIFIED_L1_INCR have not been updated yet!" );
					break;
				default:
					assert(0);	// local programming error
			}
			if(options_.get()) {
				MeritFunc_PenaltyParamUpdate_AddedStepSetOptions
					ppu_options_setter( param_update_step.get() );
				ppu_options_setter.set_options( *options_ );
			}
			merit_func_penalty_param_update_step = param_update_step;
		}

		// LineSearch_Step
		algo_step_ptr_t    line_search_full_step_step = mmp::null;
		{
			line_search_full_step_step = mmp::rcp(new LineSearchFullStep_Step(bounds_tester));
		}

		// LineSearch_Step
		algo_step_ptr_t    line_search_step = mmp::null;
		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<DirectLineSearchArmQuad_Strategy>
				direct_line_search = mmp::rcp(new  DirectLineSearchArmQuad_Strategy());
			if(options_.get()) {
				ConstrainedOptimizationPack::DirectLineSearchArmQuad_StrategySetOptions
					ls_options_setter( direct_line_search.get(), "DirectLineSearchArmQuadSQPStep" );
				ls_options_setter.set_options( *options_ );
			}
			switch( cov_.line_search_method_ ) {
				case LINE_SEARCH_DIRECT: {
					line_search_step = mmp::rcp(new LineSearchDirect_Step(direct_line_search));
					break;
				}
				case LINE_SEARCH_2ND_ORDER_CORRECT: {
					THROW_EXCEPTION(
						true, std::logic_error
						,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
						"The line_search_method option of 2ND_ORDER_CORRECT has not been updated yet!" );
					break;
				}
				case LINE_SEARCH_WATCHDOG: {
					THROW_EXCEPTION(
						true, std::logic_error
						,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
						"The line_search_method option of WATCHDOG has not been updated yet!" );
					break;
				}
			    case LINE_SEARCH_FILTER: {
					mmp::ref_count_ptr<LineSearchFilter_Step> 
						line_search_filter_step = mmp::rcp(
						  new LineSearchFilter_Step(algo_cntr->get_nlp())
						  );

					if(options_.get()) 
						{
						LineSearchFilter_StepSetOptions options_setter(line_search_filter_step.get());
						options_setter.set_options(*options_);
						}

					line_search_step = line_search_filter_step;					
				    break;
				}
			}
		}

		// LineSearchFailure
		if( new_decomp_selection_strategy.get() ) {
			line_search_step = mmp::rcp(
				new LineSearchFailureNewDecompositionSelection_Step(
					line_search_step
					,new_decomp_selection_strategy
					)
				);
		}
	
		//
		// Create the algorithm depending on the type of NLP we are trying to solve.
		//

		if( m == 0 && mI == 0 ) {

			if( nb == 0 ) {
				//
				// Unconstrained NLP (m == 0, mI == 0, num_bounded_x == 0)
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for an unconstrained "
						<< "NLP (m == 0, mI == 0, num_bounded_x == 0) ...\n";
			}
			else {
				//
				// Simple bound constrained NLP (m == 0, mI == 0, num_bounded_x > 0)
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a simple bound constrained "
						<< "NLP (m == 0, mI == 0, num_bounded_x > 0) ...\n";
			}

			int step_num       = 0;
			int assoc_step_num = 0;
	
			// EvalNewPoint
			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );
			
			// ReducedGradient
			algo->insert_step( ++step_num, ReducedGradient_name, reduced_gradient_step );

			if( nb == 0 ) {
				
				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );

			}

			// ReducedHessian
			algo->insert_step( ++step_num, ReducedHessian_name, reduced_hessian_step );

			// (.-1) Initialize reduced Hessian
			if(init_red_hess_step.get()) {
				algo->insert_assoc_step(
					step_num, IterationPack::PRE_STEP, 1
					,"InitFiniteDiffReducedHessian"
					,init_red_hess_step
					);
			}
			
			// NullSpaceStep
			algo->insert_step( ++step_num, NullSpaceStep_name, null_space_step_step );
			if( mI > 0 || nb > 0 ) {
				// SetDBoundsStd
				algo->insert_assoc_step(
					step_num
					,IterationPack::PRE_STEP
					,1
					,"SetDBoundsStd"
					,set_d_bounds_step
				  );
			}

			// CalcDFromZPZ
			algo->insert_step( ++step_num, "CalcDFromZpz", mmp::rcp(new CalcDFromZPZ_Step()) );

			if( nb > 0 ) {
				
				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );

			}

			// LineSearch
			if( cov_.line_search_method_ == LINE_SEARCH_NONE ) {
				algo->insert_step( ++step_num, LineSearch_name, line_search_full_step_step );
			}
			else {
				// Main line search step
				algo->insert_step( ++step_num, LineSearch_name, line_search_step );
				// Insert presteps
				Algorithm::poss_type
					pre_step_i = 0;
				// (.-?) LineSearchFullStep
				algo->insert_assoc_step(
					step_num
					,IterationPack::PRE_STEP
					,++pre_step_i
					,"LineSearchFullStep"
					,line_search_full_step_step
					);
				// (.-?) MeritFunc_DummyUpdate
				algo->insert_assoc_step(
					step_num
					,IterationPack::PRE_STEP
					,++pre_step_i
					,"MeritFunc_DummyUpdate"
					,mmp::rcp(new MeritFunc_DummyUpdate_Step())
					);
			}
			
		}
		else if( n == m ) {
			//
			// System of Nonlinear equations (n == m)
			//
			if(trase_out)
				*trase_out 
					<< "\nConfiguring an algorithm for a system of nonlinear equations "
					<< "NLP (n == m) ...\n";
			THROW_EXCEPTION(
				n == m, std::logic_error
				,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
				"Nonlinear equation (NLE) problems are not supported yet!" );
			assert(0); // ToDo: add the step objects for this algorithm
		}
		else if ( m > 0 || mI > 0 || nb > 0 ) {
			//
			// General nonlinear NLP ( m > 0 || mI > 0 )
			//
			if(  mI == 0 && nb == 0 ) {
				//
				// Nonlinear equality constrained NLP ( m > 0 && mI == 0 && num_bounded_x == 0 )
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a nonlinear equality constrained "
						<< "NLP ( m > 0 && mI == 0 && num_bounded_x == 0) ...\n";
			}
			else {
				//
				// Nonlinear inequality constrained NLP ( mI > 0 || num_bounded_x > 0 )
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a nonlinear generally constrained "
						<< "NLP ( mI > 0 || num_bounded_x > 0 ) ...\n";
			}

			int step_num       = 0;
			int assoc_step_num = 0;
	
			// EvalNewPoint
			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );
			if( check_descent_range_space_step_step.get() && tailored_approach && algo->algo_cntr().check_results() )
			{
				algo->insert_assoc_step(
					step_num
					,IterationPack::POST_STEP
					,1
					,"CheckDescentRangeSpaceStep"
					,check_descent_range_space_step_step
					);
			}
			
			// RangeSpaceStep
			if( !tailored_approach ) {
				algo->insert_step( ++step_num, RangeSpaceStep_name, range_space_step_step );
				assoc_step_num = 0;
				if( check_decomp_from_py_step.get() )
					algo->insert_assoc_step(
						step_num
						,IterationPack::POST_STEP
						,++assoc_step_num
						,"CheckDecompositionFromPy"
						,check_decomp_from_py_step
						);
				if( check_decomp_from_Rpy_step.get() )
					algo->insert_assoc_step(
						step_num
						,IterationPack::POST_STEP
						,++assoc_step_num
						,"CheckDecompositionFromRPy"
						,check_decomp_from_Rpy_step
						);
				if( check_descent_range_space_step_step.get() )
					algo->insert_assoc_step(
						step_num
						,IterationPack::POST_STEP
						,++assoc_step_num
						,"CheckDescentRangeSpaceStep"
						,check_descent_range_space_step_step
						);
			}

			// ReducedGradient
			if( !tailored_approach ) {
				algo->insert_step( ++step_num, ReducedGradient_name, reduced_gradient_step );
			}

			if( mI == 0 && nb == 0 ) {

				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CalcLagrangeMultDecomposed
				// Compute these here so that in case we converge we can report them
				if( !tailored_approach ) {
					// ToDo: Insert this step
				}

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
			}

			// ReducedHessian
			algo->insert_step( ++step_num, ReducedHessian_name, reduced_hessian_step );

			// (.-1) Initialize reduced Hessian
			if(init_red_hess_step.get()) {
				algo->insert_assoc_step(
					step_num, IterationPack::PRE_STEP, 1
					,"InitFiniteDiffReducedHessian"
					,init_red_hess_step
					);
			}

			// (.-1) CheckSkipBFGSUpdate
			algo->insert_assoc_step(
				step_num
				,IterationPack::PRE_STEP
				,1
				,CheckSkipBFGSUpdate_name
				,check_skip_bfgs_update_step
			  );

			// NullSpaceStep
			algo->insert_step( ++step_num, NullSpaceStep_name, null_space_step_step );
			if( mI > 0 || nb > 0 ) {
				// SetDBoundsStd
				algo->insert_assoc_step(
					step_num
					,IterationPack::PRE_STEP
					,1
					,"SetDBoundsStd"
					,set_d_bounds_step
				  );
			}

			// CalcDFromYPYZPZ
			algo->insert_step( ++step_num, CalcDFromYPYZPZ_name, calc_d_from_Ypy_Zpy_step );
			
			if( mI > 0 || nb > 0 ) {

				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CalcLagrangeMultDecomposed
				// Compute these here so that in case we converge we can report them
				if( !tailored_approach ) {
					// ToDo: Insert this step
				}

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
			}

			// LineSearch
			if( cov_.line_search_method_ == LINE_SEARCH_NONE ) {
				algo->insert_step( ++step_num, LineSearch_name, line_search_full_step_step );
			}
			else {
				// Main line search step
				algo->insert_step( ++step_num, LineSearch_name, line_search_step );
				// Insert presteps
				Algorithm::poss_type
					pre_step_i = 0;
				// (.-?) LineSearchFullStep
				algo->insert_assoc_step(
					step_num
					,IterationPack::PRE_STEP
					,++pre_step_i
					,"LineSearchFullStep"
					,line_search_full_step_step
					);
				// (.-?) MeritFunc_PenaltyPramUpdate
				if(merit_func_penalty_param_update_step.get()) {
				  algo->insert_assoc_step(
					  step_num
					  ,IterationPack::PRE_STEP
					  ,++pre_step_i
					  ,"MeritFunc_PenaltyParamUpdate"
					  ,merit_func_penalty_param_update_step
					  );
				}
			}

		}
		else {
			assert(0); // Error, this should not ever be called!
		}
	}
	
}

void rSQPAlgo_ConfigMamaJama::init_algo(rSQPAlgoInterface* _algo)
{
	using DynamicCastHelperPack::dyn_cast;
	namespace mmp = MemMngPack;

	THROW_EXCEPTION(
		_algo == NULL, std::invalid_argument
		,"rSQPAlgo_ConfigMamaJama::init_algo(_algo) : Error, "
		"_algo can not be NULL" );

	rSQPAlgo             &algo    = dyn_cast<rSQPAlgo>(*_algo);
	rSQPState	         &state   = algo.rsqp_state();
	NLP			         &nlp     = algo.nlp();
	NLPVarReductPerm     *nlp_vrp = dynamic_cast<NLPVarReductPerm*>(&nlp);
	NLPFirstOrderDirect  *nlp_fod = dynamic_cast<NLPFirstOrderDirect*>(&nlp);

	algo.max_iter( algo.algo_cntr().max_iter() );
	algo.max_run_time( algo.algo_cntr().max_run_time() );

	// Reset the iteration count to zero
	state.k(0);

	// Get organized output of vectors and matrices even if setw is not used by Step objects.
	algo.track().journal_out()
		<< std::setprecision(algo.algo_cntr().journal_print_digits())
		<< std::scientific;

	// set the first step
	algo.do_step_first(1);

	// The rest of the algorithm should initialize itself
}

// private

void rSQPAlgo_ConfigMamaJama::readin_options(
	  const OptionsFromStreamPack::OptionsFromStream     &options
	, SOptionValues                                      *ov
	, std::ostream                                       *trase_out
	)
{
	namespace	ofsp = OptionsFromStreamPack;
	using		ofsp::OptionsFromStream;
	typedef		OptionsFromStream::options_group_t		options_group_t;
	using		ofsp::StringToIntMap;
	using		ofsp::StringToBool;

	assert(ov);	// only a local class error

	// Get the options group for "rSQPAlgo_ConfigMamaJama"
	const std::string opt_grp_name = "rSQPAlgo_ConfigMamaJama";
	const OptionsFromStream::options_group_t optgrp = options.options_group( opt_grp_name );
	if( OptionsFromStream::options_group_exists( optgrp ) ) {

		// Define map for options group "MamaJama".
		const int num_opts = 11;
		enum EMamaJama {
			MAX_BASIS_COND_CHANGE_FRAC
			,EXACT_REDUCED_HESSIAN
			,QUASI_NEWTON
			,NUM_LBFGS_UPDATES_STORED
			,LBFGS_AUTO_SCALING
			,HESSIAN_INITIALIZATION
			,QP_SOLVER
			,REINIT_HESSIAN_ON_QP_FAIL
			,LINE_SEARCH_METHOD
			,MERIT_FUNCTION_TYPE
			,L1_PENALTY_PARAM_UPDATE
		};
		const char* SMamaJama[num_opts]	= {
			"max_basis_cond_change_frac"
			,"exact_reduced_hessian"
			,"quasi_newton"
			,"num_lbfgs_updates_stored"
			,"lbfgs_auto_scaling"
			,"hessian_initialization"
			,"qp_solver"
			,"reinit_hessian_on_qp_fail"
			,"line_search_method"
			,"merit_function_type"
			,"l1_penalty_parameter_update"
		};
		StringToIntMap	mama_jama_map(	opt_grp_name, num_opts, SMamaJama );

		options_group_t::const_iterator itr = optgrp.begin();
		for( ; itr != optgrp.end(); ++itr ) {
			switch( (EMamaJama)mama_jama_map( ofsp::option_name(itr) ) ) {
				case MAX_BASIS_COND_CHANGE_FRAC:
					ov->max_basis_cond_change_frac_ = ::atof( ofsp::option_value(itr).c_str() );
					break;
				case EXACT_REDUCED_HESSIAN:
					ov->exact_reduced_hessian_ = StringToBool( "exact_reduced_hessian", ofsp::option_value(itr).c_str() );
					break;
				case QUASI_NEWTON:
				{
					const std::string &opt_val = ofsp::option_value(itr);
					if( opt_val == "AUTO" )
						ov->quasi_newton_ = QN_AUTO;
					else if( opt_val == "BFGS" )
						ov->quasi_newton_ = QN_BFGS;
					else if( opt_val == "PBFGS" )
						ov->quasi_newton_ = QN_PBFGS;
					else if( opt_val == "LBFGS" )
						ov->quasi_newton_ = QN_LBFGS;
					else if( opt_val == "LPBFGS" )
						ov->quasi_newton_ = QN_LPBFGS;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"quasi_newton\" "
							", Only options of BFGS, PBFGS"
							", LBFGS, LPBFGS and AUTO are avalible."
							);
					break;
				}
				case NUM_LBFGS_UPDATES_STORED:
					ov->num_lbfgs_updates_stored_ = ::atoi( ofsp::option_value(itr).c_str() );
					break;
				case LBFGS_AUTO_SCALING:
					ov->lbfgs_auto_scaling_
						= StringToBool( "lbfgs_auto_scaling", ofsp::option_value(itr).c_str() );
					break;
				case HESSIAN_INITIALIZATION:
				{
					const std::string &opt_val = ofsp::option_value(itr);
					if( opt_val == "IDENTITY" )
						ov->hessian_initialization_ = INIT_HESS_IDENTITY;
					else if( opt_val == "FINITE_DIFF_SCALE_IDENTITY" )
						ov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_IDENTITY;
					else if( opt_val == "FINITE_DIFF_DIAGONAL" )
						ov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL;
					else if( opt_val == "FINITE_DIFF_DIAGONAL_ABS" )
						ov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS;
					else if( opt_val == "AUTO" )
						ov->hessian_initialization_ = INIT_HESS_AUTO;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"hessian_initialization\" "
							", Only options of IDENTITY, FINITE_DIFF_SCALE_IDENTITY,"
							" FINITE_DIFF_DIAGONAL, FINITE_DIFF_DIAGONAL_ABS and AUTO"
							" are available"  );
					break;
				}
				case QP_SOLVER:
				{
					const std::string &qp_solver = ofsp::option_value(itr);
					if( qp_solver == "AUTO" ) {
						ov->qp_solver_type_ = QP_AUTO;
					} else if( qp_solver == "QPSOL" ) {
						ov->qp_solver_type_ = QP_QPSOL;
					} else if( qp_solver == "QPOPT" ) {
#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
						ov->qp_solver_type_ = QP_QPOPT;
#else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : QPOPT is not supported,"
							" must define CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT!" );
#endif
					} else if( qp_solver == "QPKWIK" ) {
						ov->qp_solver_type_ = QP_QPKWIK;
					} else if( qp_solver == "QPSCHUR" ) {
						ov->qp_solver_type_ = QP_QPSCHUR;
					} else {
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"qp_solver\" "
							"Only qp solvers QPOPT, QPSOL, QPKWIK, QPSCHUR and AUTO are avalible."	);
					}
					break;
				}
				case REINIT_HESSIAN_ON_QP_FAIL:
					ov->reinit_hessian_on_qp_fail_ = StringToBool( "reinit_hessian_on_qp_fail", ofsp::option_value(itr).c_str() );
					break;
				case LINE_SEARCH_METHOD:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "NONE" ) {
						ov->line_search_method_ = LINE_SEARCH_NONE;
					} else if( option == "DIRECT" ) {
						ov->line_search_method_ = LINE_SEARCH_DIRECT;
					} else if( option == "2ND_ORDER_CORRECT" ) {
						ov->line_search_method_ = LINE_SEARCH_2ND_ORDER_CORRECT;
					} else if( option == "WATCHDOG" ) {
						ov->line_search_method_ = LINE_SEARCH_WATCHDOG;
					} else if( option == "AUTO" ) {
						ov->line_search_method_ = LINE_SEARCH_AUTO;
					} else if( option == "FILTER" ) {
						ov->line_search_method_ = LINE_SEARCH_FILTER;
					} else {
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"line_search_method\".\n"
							"Only the options NONE, DIRECT, 2ND_ORDER_CORRECT, FILTER, WATCHDOG "
							"and AUTO are avalible." );
					}
					break;
				}
				case MERIT_FUNCTION_TYPE:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "L1" )
						ov->merit_function_type_ = MERIT_FUNC_L1;
					else if( option == "MODIFIED_L1" )
						ov->merit_function_type_ = MERIT_FUNC_MOD_L1;
					else if( option == "MODIFIED_L1_INCR" )
						ov->merit_function_type_ = MERIT_FUNC_MOD_L1_INCR;
					else if( option == "AUTO" )
						ov->merit_function_type_ = MERIT_FUNC_AUTO;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"merit_function_type\".\n"
							"Only the options L1, MODIFIED_L1, MODIFIED_L1_INCR "
							"and AUTO are avalible." );
					break;
				}
				case L1_PENALTY_PARAM_UPDATE:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "WITH_MULT" )
						ov->l1_penalty_param_update_
							= L1_PENALTY_PARAM_WITH_MULT;
					else if( option == "MULT_FREE" )
						ov->l1_penalty_param_update_
							= L1_PENALTY_PARAM_MULT_FREE;
					else if( option == "AUTO" )
						ov->l1_penalty_param_update_
							= L1_PENALTY_PARAM_AUTO;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"l1_penalty_param_update\".\n"
							"Only the options WITH_MULT, MULT_FREE and AUTO"
							"are avalible."  );
					break;
				}
				default:
					assert(0);	// this would be a local programming error only.
			}
		}
	}
	else {
		if(trase_out)
			*trase_out
				<< "\n\n*** Warning!  The options group \"rSQPAlgo_ConfigMamaJama\" was not found.\n"
				<< "Using a default set of options instead ... \n";
	}
}

//
// This is where some of the default options are set and the user is alerted to what their
// value is.
//
void rSQPAlgo_ConfigMamaJama::set_default_options( 
	const SOptionValues     &uov
	,SOptionValues          *cov
	,std::ostream           *trase_out
	)
{
	if(trase_out)
		*trase_out
			<< "\n*** Setting option defaults for options not set by the user or determined some other way ...\n";

	if( cov->max_basis_cond_change_frac_ < 0.0 &&  uov.max_basis_cond_change_frac_ < 0.0 ) {
		if(trase_out)
			*trase_out
				<< "\nmax_basis_cond_change_frac < 0 : setting max_basis_cond_change_frac = 1e+4 \n";
		cov->max_basis_cond_change_frac_ = 1e+4;
	}
	else {
		cov->max_basis_cond_change_frac_ = uov.max_basis_cond_change_frac_;
	}
	cov->exact_reduced_hessian_ = uov.exact_reduced_hessian_;
	if( cov->quasi_newton_ == QN_AUTO && uov.quasi_newton_ == QN_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nquasi_newton == AUTO: setting quasi_newton = BFGS\n";
		cov->quasi_newton_ = QN_BFGS;
	}
	else if(cov->quasi_newton_ == QN_AUTO) {
		cov->quasi_newton_ = uov.quasi_newton_;
	}
	if( cov->num_lbfgs_updates_stored_ < 0 && uov.num_lbfgs_updates_stored_ < 0 ) {
		if(trase_out)
			*trase_out
				<< "\nnum_lbfgs_updates_stored < 0 : setting num_lbfgs_updates_stored = 10\n";
		cov->num_lbfgs_updates_stored_ = 10;
	}
	else if(cov->num_lbfgs_updates_stored_ < 0) {
		cov->num_lbfgs_updates_stored_ = uov.num_lbfgs_updates_stored_;
	}
	cov->lbfgs_auto_scaling_ = uov.lbfgs_auto_scaling_;
	if( cov->hessian_initialization_ == INIT_HESS_AUTO && uov.hessian_initialization_ == INIT_HESS_AUTO ) {
		if(trase_out)
 			*trase_out
				<< "\nhessian_initialization == AUTO: setting hessian_initialization = IDENTITY\n";
		cov->hessian_initialization_ = INIT_HESS_IDENTITY;
/*
		if(trase_out)
			*trase_out
				<< "\nhessian_initialization == AUTO: setting hessian_initialization = FINITE_DIFF_DIAGONAL_ABS\n";
		cov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS;
*/
	}
	else if(cov->hessian_initialization_ == INIT_HESS_AUTO) {
		cov->hessian_initialization_ = uov.hessian_initialization_;
	}
	if( cov->qp_solver_type_ == QP_AUTO && uov.qp_solver_type_ == QP_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nqp_solver_type == AUTO: setting qp_solver_type = QPKWIK\n";
		cov->qp_solver_type_ = QP_QPKWIK;
/*
		if(trase_out)
			*trase_out
				<< "\nqp_solver_type == AUTO: setting qp_solver_type = QPSCHUR\n";
		cov->qp_solver_type_ = QP_QPSCHUR;
*/
	}
	else if(cov->qp_solver_type_ == QP_AUTO) {
		cov->qp_solver_type_ = uov.qp_solver_type_;
	}
	cov->reinit_hessian_on_qp_fail_ = uov.reinit_hessian_on_qp_fail_;
	if( cov->line_search_method_ == LINE_SEARCH_AUTO && uov.line_search_method_ == LINE_SEARCH_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nline_search_method == AUTO: setting line_search_method = 2ND_ORDER_CORRECT\n";
		cov->line_search_method_ = LINE_SEARCH_2ND_ORDER_CORRECT;
	}
	else if(cov->line_search_method_ == LINE_SEARCH_AUTO) {
		cov->line_search_method_ = uov.line_search_method_;
	}
	if( cov->merit_function_type_ == MERIT_FUNC_AUTO && uov.merit_function_type_ == MERIT_FUNC_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nmerit_function_type == AUTO: setting merit_function_type = MODIFIED_L1_INCR\n";
		cov->merit_function_type_ = MERIT_FUNC_MOD_L1_INCR;
	}
	else if(cov->merit_function_type_ == MERIT_FUNC_AUTO) {
		cov->merit_function_type_ = uov.merit_function_type_;
	}
	if( cov->l1_penalty_param_update_ == L1_PENALTY_PARAM_AUTO && uov.l1_penalty_param_update_ == L1_PENALTY_PARAM_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nl1_penalty_param_update == AUTO: setting l1_penalty_param_update = MULT_FREE\n";
		cov->l1_penalty_param_update_ = L1_PENALTY_PARAM_MULT_FREE;
	}
	else if(cov->l1_penalty_param_update_ == L1_PENALTY_PARAM_AUTO) {
		cov->l1_penalty_param_update_ = uov.l1_penalty_param_update_;
	}
	if( cov->full_steps_after_k_ < 0 && uov.full_steps_after_k_ < 0 ) {
		if(trase_out)
			*trase_out
				<< "\nfull_steps_after_k < 0 : the line search will never be turned off after so many iterations\n";
	}
	else {
		cov->full_steps_after_k_ = uov.full_steps_after_k_;
	}
	if(trase_out)
		*trase_out
			<< "\n*** End setting default options\n";
}

}	// end namespace ReducedSpaceSQPPack
