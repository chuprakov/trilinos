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

#ifdef _INTEL_CXX
// disable Intel C++ 5.0 warnings about debugger limitations
#pragma warning(disable : 985)
#endif

#include <assert.h>

#include <sstream>
#include <typeinfo>
#include <iostream>

#include "Misc/include/debug.h"

#include "rSQPAlgo_ConfigMamaJama.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgoContainer.h"
//#include "ReducedSpaceSQPPack/include/rSQPStateContinuousStorage.h"
//#include "ReducedSpaceSQPPack/include/rSQPStateContinuousStorageMatrixWithOpCreatorAggr.h"
//#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
//#include "SparseLinAlgPack/include/COOMatrixWithPartitionedViewSubclass.h"			// HL, Gc, Hcj
#include "ConstrainedOptimizationPack/include/MatrixIdentConcatStd.h"                   // Y, Z
//#include "SparseLinAlgPack/include/COOMatrixPartitionViewSubclass.h"				// U
//#include "SparseLinAlgPack/include/GenMatrixSubclass.h"								// V
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefCholFactor.h"          // rHL 
//#include "ConstrainedOptimizationPack/include/MatrixSymPosDefInvCholFactor.h"		// .
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefLBFGS.h"				// .
//#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasicInitDiagonal.h"// | rHL (super basics)
//#include "SparseLinAlgPack/include/MatrixSymDiagonalStd.h"                          // |

#include "ConstrainedOptimizationPack/include/VariableBoundsTesterSetOptions.h"

#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "NLPInterfacePack/test/NLPFirstOrderDirectTester.h"
#include "NLPInterfacePack/test/NLPFirstOrderDirectTesterSetOptions.h"

#include "NLPInterfacePack/include/CalcFiniteDiffProd.h"
#include "NLPInterfacePack/include/CalcFiniteDiffProdSetOptions.h"
#include "NLPInterfacePack/test/NLPFirstDerivativesTester.h"
#include "NLPInterfacePack/test/NLPFirstDerivativesTesterSetOptions.h"

#include "NLPInterfacePack/include/NLPSecondOrderInfo.h"
#include "NLPInterfacePack/include/NLPReduced.h"
//#include "SparseSolverPack/include/COOBasisSystem.h"
#ifdef SPARSE_SOLVER_PACK_USE_MA28
//#include "SparseSolverPack/include/MA28SparseCOOSolverCreator.h"
#endif
#ifdef SPARSE_SOLVER_PACK_USE_MA48
//#include "SparseSolverPack/include/MA48SparseCOOSolverCreator.h"
#endif

// line search
#include "ConstrainedOptimizationPack/include/DirectLineSearchArmQuad_Strategy.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearchArmQuad_StrategySetOptions.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPL1.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPModL1.h"

#include "ConstrainedOptimizationPack/include/DecompositionSystemTester.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemTesterSetOptions.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinate.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemOrthogonal.h"

#include "AbstractLinAlgPack/include/BasisSystemTester.h"
#include "AbstractLinAlgPack/include/BasisSystemTesterSetOptions.h"

//#include "ConstrainedOptimizationPack/include/QPSolverRelaxedTester.h"
//#include "ConstrainedOptimizationPack/include/QPSolverRelaxedTesterSetOptions.h"
//#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPSchur.h"
//#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPSchurSetOptions.h"
//#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianFull.h"
//#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianSuperBasic.h"
//#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPKWIK.h"
#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
//#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPOPT.h"
#endif

#include "ReducedSpaceSQPPack/include/std/rSQPAlgorithmStepNames.h"

#include "ReducedSpaceSQPPack/include/std/EvalNewPointStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachCoordinate_Step.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachOrthogonal_Step.h"
#include "ReducedSpaceSQPPack/include/std/ReducedGradientStd_Step.h"
//#include "ReducedSpaceSQPPack/include/std/InitFinDiffReducedHessian_Step.h"
//#include "ReducedSpaceSQPPack/include/std/InitFinDiffReducedHessian_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSFull_Strategy.h"
//#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSProjected_Strategy.h"
//#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.h"
//#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateLPBFGS_Strategy.h"
//#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.h"
#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_StrategySetOptions.h"
#include "ReducedSpaceSQPPack/include/std/RangeSpaceStepStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/CheckDescentRangeSpaceStep_Step.h"
//#include "ReducedSpaceSQPPack/include/std/CheckBasisFromCPy_Step.h"
//#include "ReducedSpaceSQPPack/include/std/CheckBasisFromPy_Step.h"
#include "ReducedSpaceSQPPack/include/std/NullSpaceStepWithoutBounds_Step.h"
//#include "ReducedSpaceSQPPack/include/std/SetDBoundsStd_AddedStep.h"
//#include "ReducedSpaceSQPPack/include/std/QPFailureReinitReducedHessian_Step.h"
//#include "ReducedSpaceSQPPack/include/std/IndepDirecWithBoundsStd_Step.h"
//#include "ReducedSpaceSQPPack/include/std/IndepDirecWithBoundsStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/CalcDFromYPYZPZ_Step.h"
//#include "ReducedSpaceSQPPack/include/std/LineSearchFailureNewBasisSelection_Step.h"
//#include "ReducedSpaceSQPPack/include/std/NewBasisSelectionStd_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchFullStep_Step.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchDirect_Step.h"
//#include "ReducedSpaceSQPPack/include/std/LineSearch2ndOrderCorrect_Step.h"
//#include "ReducedSpaceSQPPack/include/std/LineSearch2ndOrderCorrect_StepSetOptions.h"
//#include "ReducedSpaceSQPPack/include/std/FeasibilityStepReducedStd_Strategy.h"
//#include "ReducedSpaceSQPPack/include/std/FeasibilityStepReducedStd_StrategySetOptions.h"
//#include "ReducedSpaceSQPPack/include/std/QuasiRangeSpaceStepStd_Strategy.h"
//#include "ReducedSpaceSQPPack/include/std/QuasiRangeSpaceStepTailoredApproach_Strategy.h"
//#include "ReducedSpaceSQPPack/include/std/LineSearchWatchDog_Step.h"
//#include "ReducedSpaceSQPPack/include/std/LineSearchWatchDog_StepSetOptions.h"
//#include "ReducedSpaceSQPPack/include/std/LineSearchFullStepAfterKIter_Step.h"
//#include "ReducedSpaceSQPPack/include/std/CalcLambdaIndepStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CalcReducedGradLagrangianStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CheckConvergenceStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CheckConvergenceStd_AddedStepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/CheckSkipBFGSUpdateStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.h"
//#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateWithMult_AddedStep.h"
//#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.h"
//#include "ReducedSpaceSQPPack/include/std/MeritFunc_ModifiedL1LargerSteps_AddedStep.h"
//#include "ReducedSpaceSQPPack/include/std/MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.h"
//#include "ReducedSpaceSQPPack/include/std/ActSetStats_AddedStep.h"
//#include "ReducedSpaceSQPPack/include/std/NumFixedDepIndep_AddedStep.h"

//#include "ReducedSpaceSQPPack/include/std/act_set_stats.h"
//#include "ReducedSpaceSQPPack/include/std/qp_solver_stats.h"
//#include "ReducedSpaceSQPPack/include/std/quasi_newton_stats.h"

//#include "SparseLinAlgPack/include/sparse_bounds.h"

//#include "LinAlgPack/include/PermVecMat.h"

// Misc utilities
#include "AbstractFactoryStd.h"
#include "dynamic_cast_verbose.h"
#include "ReleaseResource_ref_count_ptr.h"
#include "ThrowException.h"

// Stuff to read in options
#include "Misc/include/StringToIntMap.h"
#include "Misc/include/StringToBool.h"

// Stuff for exact reduced hessian
//#include "ReducedSpaceSQPPack/include/std/ReducedHessianExactStd_Step.h"
//#include "ReducedSpaceSQPPack/include/std/CrossTermExactStd_Step.h"
//#include "ReducedSpaceSQPPack/include/std/DampenCrossTermStd_Step.h"

namespace {
	const double INF_BASIS_COND_CHANGE_FRAC      = 1e+20;
	const int DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE = 200;
}

namespace ReducedSpaceSQPPack {

//
// Here is where we define the default values for the algorithm.  These
// should agree with what are in the rSQPpp.opt.rSQPAlgo_ConfigMamaJama file.
//
rSQPAlgo_ConfigMamaJama::SOptionValues::SOptionValues()
	:direct_linear_solver_type_(LA_AUTO)
	,null_space_matrix_type_(NULL_SPACE_MATRIX_AUTO)
	,range_space_matrix_type_(RANGE_SPACE_MATRIX_AUTO)
	,max_basis_cond_change_frac_(-1.0)
	,exact_reduced_hessian_(false)
	,quasi_newton_(QN_AUTO)
	,max_dof_quasi_newton_dense_(-1)
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

rSQPAlgo_ConfigMamaJama::rSQPAlgo_ConfigMamaJama(
	const basis_sys_ptr_t                     &basis_sys
	,const var_reduct_orthog_strategy_ptr_t   &var_reduct_orthog_strategy
	)
{
	this->initialize(basis_sys,var_reduct_orthog_strategy);
}

void rSQPAlgo_ConfigMamaJama::initialize(
	const basis_sys_ptr_t                     &basis_sys
	,const var_reduct_orthog_strategy_ptr_t   &var_reduct_orthog_strategy
	)
{
	basis_sys_                  = basis_sys;
	var_reduct_orthog_strategy_ = var_reduct_orthog_strategy;
}

rSQPAlgo_ConfigMamaJama::~rSQPAlgo_ConfigMamaJama()
{
	// No need to really do anything!
}

// overridden from rSQPAlgo_Config

void rSQPAlgo_ConfigMamaJama::set_options( const options_ptr_t& options )
{
	options_ = options;
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
	namespace afp = AbstractFactoryPack;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
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

	typedef rcp::ref_count_ptr<rSQPAlgo>	algo_ptr_t;
	algo_ptr_t algo = rcp::rcp(new rSQPAlgo);
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

	if(trase_out)
		*trase_out << "\n*** Probing the NLP object for supported interfaces ...\n";

	// Get the dimensions of the NLP
	NLP &nlp = algo->nlp();
	if(!nlp.is_initialized()) nlp.initialize();
	const size_type
		n   = algo->nlp().n(),
		m   = algo->nlp().m(),
		mI  = algo->nlp().mI(),
		r   = m, // ToDo: Compute this for real!
		dof = n - r,
		nb  = nlp.num_bounded_x();

	// Determine which NLP interface is supported
	NLPFirstOrderInfo    *nlp_foi = dynamic_cast<NLPFirstOrderInfo*>(   algo->get_nlp() );	
	NLPSecondOrderInfo   *nlp_soi = dynamic_cast<NLPSecondOrderInfo*>(  algo->get_nlp() );	
	NLPFirstOrderDirect  *nlp_fod = dynamic_cast<NLPFirstOrderDirect*>( algo->get_nlp() );
	bool tailored_approach = false;
	if( nlp_foi ) {
		if(trase_out)
			*trase_out << "\nDetected that NLP object supports the NLPFirstOrderInfo interface!\n";
		tailored_approach = false;
	}
	else {
		if( nlp_fod ) {
			if(trase_out)
				*trase_out << "\nDetected that NLP object supports the NLPFirstOrderDirect interface!\n";
			tailored_approach = true;
		}
		else {
			THROW_EXCEPTION(
				true, std::logic_error
				,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
				"the NLP object of type \'" << typeid(algo->nlp()).name() <<
				"\' does not support the NLPFirstOrderDirect or "
				"NLPFirstOrderInfo interfaces!" );
		}
	}
	if( nlp_soi ) {
		if(trase_out)
			*trase_out << "\nDetected that NLP object also supports the NLPSecondOrderInfo interface!\n";
	}

	// //////////////////////////////////////////////////////
	// C.1.  Sort out the options

	if(trase_out)
		*trase_out
			<< "\n*** Sorting out some of the options given input options ...\n";

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
		cov_.merit_function_type_		= MERIT_FUNC_L1;
		cov_.l1_penalty_param_update_	= L1_PENALTY_PARAM_MULT_FREE;
		cov_.null_space_matrix_type_    = NULL_SPACE_MATRIX_EXPLICIT;
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

	// Set default
	if( uov_.max_dof_quasi_newton_dense_ < 0 )
		cov_.max_dof_quasi_newton_dense_ = DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE;
	else
		cov_.max_dof_quasi_newton_dense_ = uov_.max_dof_quasi_newton_dense_;

	// Decide what type of quasi newton update to use
	switch( uov_.quasi_newton_ ) {
		case QN_AUTO: {
			if(trase_out)
				*trase_out
					<< "\nquasi_newton == AUTO:"
					<< "\nnlp.num_bounded_x() == " << nlp.num_bounded_x() << ":\n";
			if( n - r > cov_.max_dof_quasi_newton_dense_ ) {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " > max_dof_quasi_newton_dense = "
						<< cov_.max_dof_quasi_newton_dense_ <<  ":\n"
						<< "setting quasi_newton == LBFGS\n";
				cov_.quasi_newton_ = QN_LBFGS;
			}
			else {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " <= max_dof_quasi_newton_dense = "
						<< cov_.max_dof_quasi_newton_dense_ << ":\n"
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

	// Decide what type of range space matrix to use
	if( uov_.range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO ) {
		const bool use_orth = dof*dof*r	<= cov_.max_dof_quasi_newton_dense_*cov_.max_dof_quasi_newton_dense_;
		if(trase_out)
			*trase_out
				<< "\nrange_space_matrix == AUTO:"
				<< "\n(n-r)^2*r = (" << dof << ")^2 * " << r << " = " << (dof*dof*r)
				<< ( use_orth ? " <= " : " > " ) << "max_dof_quasi_newton_dense^2 = ("
				<< cov_.max_dof_quasi_newton_dense_ << ")^2 = "
				<< cov_.max_dof_quasi_newton_dense_*cov_.max_dof_quasi_newton_dense_
				<< ( use_orth
					 ? "\nsetting range_space_matrix = ORTHOGONAL\n"
					 : "\nsetting range_space_matrix = COORDINATE\n" );
		cov_.range_space_matrix_type_ =
			( use_orth
			  ? RANGE_SPACE_MATRIX_ORTHOGONAL
			  : RANGE_SPACE_MATRIX_COORDINATE );
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

	// /////////////////////////////////////////////////////
	// C.1. Create the decomposition system object

	typedef ref_count_ptr<DecompositionSystem> decomp_sys_ptr_t;
	decomp_sys_ptr_t decomp_sys = rcp::null;
	if(!tailored_approach) {
		// Set the default basis system if one is not set
		THROW_EXCEPTION(
			basis_sys_.get() == NULL, std::logic_error
			,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
			"There is no default basis system yet!" );
		// Create the testing object for the basis system and set it up.
		ref_count_ptr<BasisSystemTester>
			basis_sys_tester = rcp::rcp(new BasisSystemTester());
		if(options_.get()) {
			BasisSystemTesterSetOptions
				opt_setter(basis_sys_tester.get());
			opt_setter.set_options(*options_);
		}
		// create the step object
		switch( cov_.range_space_matrix_type_ ) {
			case RANGE_SPACE_MATRIX_COORDINATE:
				decomp_sys
					= rcp::rcp(new DecompositionSystemCoordinate(
						nlp.space_x()
						,nlp.space_c()
						,nlp.space_h()
						,basis_sys_
						,basis_sys_tester
						) );
				break;
			case RANGE_SPACE_MATRIX_ORTHOGONAL: {
				// Create a default VarReductOrthog_Strategy object for serial applications!
				if( var_reduct_orthog_strategy_.get() == NULL )
				{
					assert(0); // ToDo: Implement
				}
				decomp_sys
					= rcp::rcp(new DecompositionSystemOrthogonal(
						nlp.space_x()
						,nlp.space_c()
						,nlp.space_h()
						,basis_sys_
						,basis_sys_tester
						,var_reduct_orthog_strategy_
						) );
				break;
			}
			default:
				assert(0);	// only a local error
		}
	}

	// /////////////////////////////////////////////////////
	// C.2. Create and set the state object

	if(trase_out)
		*trase_out
			<< "\n*** Creating the state object an setting up iteration quantity objects ...\n";

	{
		//
		// Create the state object with the vector spaces
		//

		typedef ref_count_ptr<rSQPState>   state_ptr_t;
		state_ptr_t
			state = rcp::rcp(
				new rSQPState(
					decomp_sys
					,nlp.space_x()
					,nlp.space_c()
					,nlp.space_h()
					,(tailored_approach
					  ? ( nlp_fod->var_dep().size() 
					    ? nlp.space_x()->sub_space(nlp_fod->var_dep())->clone()
					    : rcp::null )
                      : decomp_sys->space_range() )
					,(tailored_approach
					   ?( nlp_fod->var_indep().size()
					     ? nlp.space_x()->sub_space(nlp_fod->var_indep())->clone()
					     : rcp::null )
					   : decomp_sys->space_null() )
					)
				);

		//
		// Set the iteration quantities for the NLP matrix objects
		//

		THROW_EXCEPTION( // ToDo: Remove this and support mI > 0 in the futrue!
			mI, std::logic_error
			,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
			"general inequaity constraints are not supported yet!" );
		if( tailored_approach ) {
			// NLPFirstOrderDirect
			assert( mI == 0 && nlp_fod->con_undecomp().size() == 0 );
			// ToDo: Add the necessary iteration quantities when mI > 0 and
			// con_undecomp().size() > 0 are supported!
		}
		else {
			// NLPFirstOrderInfo
			if(m)
				state->set_iter_quant(
					Gc_name
					,rcp::rcp(
						new IterQuantityAccessContiguous<MatrixWithOp>(
							1
							,Gc_name
							,nlp_foi->space_Gc()
							)
						)
					);
			if(mI)
				state->set_iter_quant(
					Gh_name
					,rcp::rcp(
						new IterQuantityAccessContiguous<MatrixWithOp>(
							1
							,Gh_name
							,nlp_foi->space_Gh()
							)
						)
					);
			if(nlp_soi)
				state->set_iter_quant(
					HL_name
					,rcp::rcp(
						new IterQuantityAccessContiguous<MatrixSymWithOp>(
							1
							,HL_name
							,nlp_soi->space_HL()
							)
						)
					);
		}

		//
		// Set the algorithm specific matrix objects
		//
		
		// Add range/null decomposition matrices

		if(tailored_approach) {
			// Z
			state->set_iter_quant(
				Z_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<MatrixWithOp>(
						1
						,Z_name
						,rcp::rcp(new afp::AbstractFactoryStd<MatrixWithOp,MatrixIdentConcatStd>)
						)
					)
				);
			// Y
			state->set_iter_quant(
				Y_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<MatrixWithOp>(
						1
						,Y_name
						,rcp::rcp(new afp::AbstractFactoryStd<MatrixWithOp,MatrixIdentConcatStd>)
						)
					)
				);
			// ToDo: Add matrix iq object for Uz
			// ToDo: Add matrix iq object for Uy
			// ToDo: Add matrix iq object for Vz
			// ToDo: Add matrix iq object for Vy
		}
		else {
			// Z
			state->set_iter_quant(
				Z_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<MatrixWithOp>(
						1
						,Z_name
						,decomp_sys->factory_Z()
						)
					)
				);
			// Y
			state->set_iter_quant(
				Y_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<MatrixWithOp>(
						1
						,Y_name
						,decomp_sys->factory_Y()
						)
					)
				);
			// R
			state->set_iter_quant(
				R_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<MatrixWithOpNonsingular>(
						1
						,R_name
						,decomp_sys->factory_R()
						)
					)
				);
			// ToDo: Add matrix iq object for Uz
			// ToDo: Add matrix iq object for Uy
			// ToDo: Add matrix iq object for Vz
			// ToDo: Add matrix iq object for Vy
		}

		// Add reduced Hessian

		if( !cov_.exact_reduced_hessian_ ) {
			ref_count_ptr<afp::AbstractFactory<MatrixSymWithOp> >
				abstract_factory_rHL = rcp::null;
			// Only maintain the orginal matrix if we have inequality constraints and therefore will be
			// needing a QP solver (which may be QPSchur which needs an accurate original matrix for
			// iterative refinment).
			const bool
				maintain_original = ( nb || mI );
			// Maintain the factor if a QP solver is needed, QPSchur is being used, or we are checking
			// results
			const bool
				maintain_inverse = ( (!nb && !mI && m==r) || cov_.qp_solver_type_==QP_QPSCHUR
									 || algo->algo_cntr().check_results() );
			switch( cov_.quasi_newton_ ) {
				case QN_BFGS:
					abstract_factory_rHL = rcp::rcp(
						new afp::AbstractFactoryStd<MatrixSymWithOp,MatrixSymPosDefCholFactor,MatrixSymPosDefCholFactor::PostMod>(
							MatrixSymPosDefCholFactor::PostMod(
								maintain_original      // maintain_original
								,maintain_inverse      // maintain_factor
								,true                  // allow_factor      (always!)
								)
							)
						);
					break;
				case QN_LBFGS:
					abstract_factory_rHL = rcp::rcp(
						new afp::AbstractFactoryStd<MatrixSymWithOp,MatrixSymPosDefLBFGS,MatrixSymPosDefLBFGS::PostMod>(
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
				,rcp::rcp(
					new IterQuantityAccessContiguous<MatrixSymWithOp>(
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

		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<afp::AbstractFactory<MeritFuncNLP> >
				merit_func_factory = rcp::null;
			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1:
					merit_func_factory = rcp::rcp(
						new afp::AbstractFactoryStd<MeritFuncNLP,MeritFuncNLPL1>());
					break;
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
					merit_func_factory = rcp::rcp(
						new afp::AbstractFactoryStd<MeritFuncNLP,MeritFuncNLPModL1>());
					break;
				default:
					assert(0);	// local programming error
			}
			state->set_iter_quant(
				merit_func_nlp_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<MeritFuncNLP>(
						1
						,merit_func_nlp_name
						,merit_func_factory
						)
					)
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
		typedef IterQuantityAccessContiguous<VectorWithOpMutable>   IQ_vector_cngs;

		dyn_cast<IQ_vector_cngs>(state->x()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->f()).resize(2);
		if(m) dyn_cast<IQ_vector_cngs>(state->c()).resize(2);
		state->Gf();
		if(mI) dyn_cast<IQ_vector_cngs>(state->h()).resize(2);
		if(m && nlp_foi) state->Gc();
		if(mI) state->Gh();

		if(m) state->py();
		if(m) dyn_cast<IQ_vector_cngs>(state->Ypy()).resize(2);
		if(m) state->pz();
		if(m) dyn_cast<IQ_vector_cngs>(state->Zpz()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->d()).resize(2);

		dyn_cast<IQ_vector_cngs>(state->rGf()).resize(2);
		state->w();
		state->zeta();
		state->qp_grad();
		state->eta();

		dyn_cast<IQ_scalar_cngs>(state->alpha()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->mu()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->phi()).resize(2);

		dyn_cast<IQ_scalar_cngs>(state->opt_kkt_err()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->feas_kkt_err()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->rGL()).resize(2);
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
		
		typedef ref_count_ptr<AlgorithmStep>    algo_step_ptr_t;

		// Create the variable bounds testing object.
		typedef rcp::ref_count_ptr<VariableBoundsTester>     bounds_tester_ptr_t;
		bounds_tester_ptr_t   bounds_tester = rcp::null;
		if(nb) { // has variable bounds?
			const value_type var_bounds_warning_tol = 1e-10;
			const value_type var_bounds_error_tol   = 1e-5;
			bounds_tester = rcp::rcp(
				new VariableBoundsTester(
					var_bounds_warning_tol      // default warning tolerance
					,var_bounds_error_tol       // default error tolerance
					) );
			if(options_.get()) {
				ConstrainedOptimizationPack::VariableBoundsTesterSetOptions
					options_setter( bounds_tester.get() );
				options_setter.set_options(*options_);
				}
		}

		// Create the finite difference class
		typedef rcp::ref_count_ptr<CalcFiniteDiffProd>     calc_fd_prod_ptr_t;
		calc_fd_prod_ptr_t   calc_fd_prod = rcp::null;
		{
			calc_fd_prod = rcp::rcp(new CalcFiniteDiffProd());
			if(options_.get()) {
				ConstrainedOptimizationPack::CalcFiniteDiffProdSetOptions
					options_setter( calc_fd_prod.get() );
				options_setter.set_options(*options_);
			}
		}
		
		// EvalNewPoint_Step
		algo_step_ptr_t    eval_new_point_step = rcp::null;
		{
			// Create the step object
			if( tailored_approach ) {
				// create and setup the derivative tester
				typedef rcp::ref_count_ptr<NLPFirstOrderDirectTester>   deriv_tester_ptr_t;
				deriv_tester_ptr_t
					deriv_tester = rcp::rcp(
						new NLPFirstOrderDirectTester(
							calc_fd_prod
							,NLPFirstOrderDirectTester::FD_DIRECTIONAL    // Gf testing
							,NLPFirstOrderDirectTester::FD_DIRECTIONAL    // -Inv(C)*N testing
							) );
				if(options_.get()) {
					NLPInterfacePack::NLPFirstOrderDirectTesterSetOptions
						options_setter(deriv_tester.get());
					options_setter.set_options(*options_);
				}
				// create the step
				typedef rcp::ref_count_ptr<EvalNewPointTailoredApproach_Step>  _eval_new_point_step_ptr_t;
				_eval_new_point_step_ptr_t
					_eval_new_point_step = rcp::null;
				switch( cov_.range_space_matrix_type_ ) {
					case RANGE_SPACE_MATRIX_COORDINATE:
						_eval_new_point_step
							= rcp::rcp(new EvalNewPointTailoredApproachCoordinate_Step(deriv_tester,bounds_tester));
						break;
					case RANGE_SPACE_MATRIX_ORTHOGONAL:
						_eval_new_point_step
							= rcp::rcp(new EvalNewPointTailoredApproachOrthogonal_Step(
								var_reduct_orthog_strategy_,deriv_tester,bounds_tester) );
						break;
					default:
						assert(0);	// only a local error
				}
				if(options_.get()) {
					EvalNewPointTailoredApproach_StepSetOptions
						options_setter(_eval_new_point_step.get());
					options_setter.set_options(*options_);
				}
				eval_new_point_step = _eval_new_point_step;
			}
			else {
				// create and setup the derivative tester
				typedef rcp::ref_count_ptr<NLPFirstDerivativesTester>   deriv_tester_ptr_t;
				deriv_tester_ptr_t
					deriv_tester = rcp::rcp(
						new NLPFirstDerivativesTester(
							calc_fd_prod
							,NLPFirstDerivativesTester::FD_DIRECTIONAL
							) );
				if(options_.get()) {
					NLPInterfacePack::NLPFirstDerivativesTesterSetOptions
						options_setter(deriv_tester.get());
					options_setter.set_options(*options_);
				}
				// create and setup the decomposition system tester
				typedef rcp::ref_count_ptr<DecompositionSystemTester>   decomp_sys_tester_ptr_t;
				decomp_sys_tester_ptr_t
					decomp_sys_tester = rcp::rcp( new DecompositionSystemTester() );
				if(options_.get()) {
					DecompositionSystemTesterSetOptions
						options_setter(decomp_sys_tester.get());
					options_setter.set_options(*options_);
				}
				typedef rcp::ref_count_ptr<EvalNewPointStd_Step>  _eval_new_point_step_ptr_t;
				_eval_new_point_step_ptr_t
					_eval_new_point_step = rcp::rcp(
						new EvalNewPointStd_Step(
							deriv_tester
							,decomp_sys_tester
							,bounds_tester
							) );
				if(options_.get()) {
					EvalNewPointStd_StepSetOptions
						options_setter(_eval_new_point_step.get());
					options_setter.set_options(*options_);
				}
				eval_new_point_step = _eval_new_point_step;
			}
		}

		// RangeSpace_Step
		algo_step_ptr_t    range_space_step_step = rcp::null;
		if( !tailored_approach ) {
			range_space_step_step = rcp::rcp(new RangeSpaceStepStd_Step());
		}

		// CheckDescentRangeSpaceStep
		algo_step_ptr_t    check_descent_range_space_step_step = rcp::null;
		if( algo->algo_cntr().check_results() ) {
			check_descent_range_space_step_step = rcp::rcp(new CheckDescentRangeSpaceStep_Step(calc_fd_prod));
		}

		// ReducedGradient_Step
		algo_step_ptr_t    reduced_gradient_step = rcp::null;
		if( !tailored_approach ) {
			reduced_gradient_step = rcp::rcp(new ReducedGradientStd_Step());
		}

		// CheckSkipBFGSUpdate
		algo_step_ptr_t    check_skip_bfgs_update_step = rcp::null;
		if(!cov_.exact_reduced_hessian_) {
			ref_count_ptr<CheckSkipBFGSUpdateStd_Step>
				step = rcp::rcp(new CheckSkipBFGSUpdateStd_Step());
			if(options_.get()) {
				CheckSkipBFGSUpdateStd_StepSetOptions
					opt_setter( step.get() );
				opt_setter.set_options( *options_ );
			}
			check_skip_bfgs_update_step = step;
		}

		// ReducedHessian_Step
		algo_step_ptr_t    reduced_hessian_step = rcp::null;
		{
			// Get the strategy object that will perform the actual secant update.
			ref_count_ptr<ReducedHessianSecantUpdate_Strategy>
				secant_update_strategy = rcp::null;
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
						bfgs_strategy = rcp::rcp(new BFGSUpdate_Strategy);
					if(options_.get()) { 
						BFGSUpdate_StrategySetOptions
							opt_setter( bfgs_strategy.get() );
						opt_setter.set_options( *options_ );
					}
					switch( cov_.quasi_newton_ ) {
					    case QN_BFGS:
					    case QN_LBFGS:
						{
							secant_update_strategy = rcp::rcp(new ReducedHessianSecantUpdateBFGSFull_Strategy(bfgs_strategy));
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
			reduced_hessian_step = rcp::rcp(
				new ReducedHessianSecantUpdateStd_Step( secant_update_strategy ) );
			// Add the QuasiNewtonStats iteration quantity
			algo->state().set_iter_quant(
				quasi_newton_stats_name
				,rcp::rcp(new IterQuantityAccessContiguous<QuasiNewtonStats>(
					1
					,quasi_newton_stats_name
#ifdef _MIPS_CXX
					,rcp::ref_count_ptr<AbstractFactoryPack::AbstractFactoryStd<QuasiNewtonStats,QuasiNewtonStats> >(
						new AbstractFactoryPack::AbstractFactoryStd<QuasiNewtonStats,QuasiNewtonStats>())
#endif
					)
				));
		}

		// NullSpace_Step
		algo_step_ptr_t    null_space_step_step = rcp::null;
		{
			null_space_step_step = rcp::rcp(new NullSpaceStepWithoutBounds_Step());
		}

		// CalcDFromYPYZPZ_Step
		algo_step_ptr_t    calc_d_from_Ypy_Zpy_step = rcp::null;
		{
			calc_d_from_Ypy_Zpy_step = rcp::rcp(new CalcDFromYPYZPZ_Step());
		}

		// CalcReducedGradLagrangianStd_AddedStep
		algo_step_ptr_t    calc_reduced_grad_lagr_step = rcp::null;
		{
			calc_reduced_grad_lagr_step = rcp::rcp(
				new CalcReducedGradLagrangianStd_AddedStep() );
		}

		// CheckConvergence_Step
		algo_step_ptr_t    check_convergence_step = rcp::null;
		{
			ref_count_ptr<CheckConvergenceStd_AddedStep>
				_check_convergence_step = rcp::rcp(new CheckConvergenceStd_AddedStep());
			if(options_.get()) {
				CheckConvergenceStd_AddedStepSetOptions
					opt_setter( _check_convergence_step.get() );
				opt_setter.set_options( *options_ );
			}
			check_convergence_step = _check_convergence_step;
		}

		// MeritFuncPenaltyParamUpdate_Step
		algo_step_ptr_t    merit_func_penalty_param_update_step = rcp::null;
		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<MeritFunc_PenaltyParamUpdate_AddedStep>
				param_update_step = rcp::null;
			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1: {
					switch(cov_.l1_penalty_param_update_) {
						case L1_PENALTY_PARAM_WITH_MULT:
//							param_update_step
//								= rcp::rcp(new  MeritFunc_PenaltyParamUpdateWithMult_AddedStep());
							THROW_EXCEPTION(
								true, std::logic_error
								,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
								"The l1_penalty_parameter_update option of MULT_FREE has not been updated yet!" );
							break;
						case L1_PENALTY_PARAM_MULT_FREE:
							param_update_step
								= rcp::rcp(new  MeritFunc_PenaltyParamUpdateMultFree_AddedStep());
							break;
						default:
							assert(0);
					}
					break;
				}
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
//					param_update_step = new  MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(
//											rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
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
		algo_step_ptr_t    line_search_full_step_step = rcp::null;
		{
			line_search_full_step_step = rcp::rcp(new LineSearchFullStep_Step(bounds_tester));
		}

		// LineSearch_Step
		algo_step_ptr_t    line_search_step = rcp::null;
		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<DirectLineSearchArmQuad_Strategy>
				direct_line_search = rcp::rcp(new  DirectLineSearchArmQuad_Strategy());
			if(options_.get()) {
				ConstrainedOptimizationPack::DirectLineSearchArmQuad_StrategySetOptions
					ls_options_setter( direct_line_search.get(), "DirectLineSearchArmQuadSQPStep" );
				ls_options_setter.set_options( *options_ );
			}
			switch( cov_.line_search_method_ ) {
				case LINE_SEARCH_DIRECT: {
					line_search_step = rcp::rcp(new LineSearchDirect_Step(direct_line_search));
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
			}
		}
		
		//
		// Create the algorithm depending on the type of NLP we are
		// trying to solve.
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
				THROW_EXCEPTION(
					m == 0 && mI == 0 && nb == 0, std::logic_error
					,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
					"Unconstrained NLPs are not supported yet!" );
			}
			else {
				//
				// Simple bound constrained NLP (m == 0, mI == 0, num_bounded_x > 0)
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a simple bound constrained "
						<< "NLP (m == 0, mI == 0, num_bounded_x > 0) ...\n";
				THROW_EXCEPTION(
					m == 0 && mI == 0 && nb == 0, std::logic_error
					,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
					"Bound constrained NLPs are not supported yet!" );
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
		else if ( m > 0 && mI == 0 && nb == 0 ) {
			//
			// Nonlinear equality constrained NLP ( m > 0 && mI == 0 && num_bounded_x == 0 )
			//
			if(trase_out)
				*trase_out 
					<< "\nConfiguring an algorithm for nonlinear equality constrained "
					<< "NLP ( m > 0 && mI == 0 && num_bounded_x == 0) ...\n";

			int step_num = 0;
	
			// (1) EvalNewPoint
			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );
			if( tailored_approach && algo->algo_cntr().check_results() ) {
				algo->insert_assoc_step(
					step_num
					,GeneralIterationPack::POST_STEP
					,1
					,"CheckDescentRangeSpaceStep"
					,check_descent_range_space_step_step
					);
			}
			
			// (2) RangeSpaceStep
			if( !tailored_approach ) {
				algo->insert_step( ++step_num, RangeSpaceStep_name, range_space_step_step );
				if(algo->algo_cntr().check_results() ) {
					algo->insert_assoc_step(
						step_num
						,GeneralIterationPack::POST_STEP
						,1
						,"CheckDescentRangeSpaceStep"
						,check_descent_range_space_step_step
						);
				}
			}

			// (3) ReducedGradient
			if( !tailored_approach ) {
				algo->insert_step( ++step_num, ReducedGradient_name, reduced_gradient_step );
			}

			// (4) CalcReducedGradLagrangian
			algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

			// (5) CalcLagrangeMultDecomposed
			// Compute these here so that in case we converge we can report them
			if( !tailored_approach ) {
//				assert(0); // ToDo: Insert this step
			}

			// (6) CheckConvergence
			algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );

			// (7) ReducedHessian
			algo->insert_step( ++step_num, ReducedHessian_name, reduced_hessian_step );

			// (7.-1) CheckSkipBFGSUpdate
			algo->insert_assoc_step(
				step_num
				,GeneralIterationPack::PRE_STEP
				,1
				,CheckSkipBFGSUpdate_name
				,check_skip_bfgs_update_step
			  );

			// (8) NullSpaceStep
			algo->insert_step( ++step_num, NullSpaceStep_name, null_space_step_step );

			// (9) CalcDFromYPYZPZ
			algo->insert_step( ++step_num, CalcDFromYPYZPZ_name, calc_d_from_Ypy_Zpy_step );
			// ToDo: Add ths step!
			
			// (10) LineSearch
			if( cov_.line_search_method_ == LINE_SEARCH_NONE ) {
				algo->insert_step( ++step_num, LineSearch_name, line_search_full_step_step );
			}
			else {
				// (10) Main line search step
				algo->insert_step( ++step_num, LineSearch_name, line_search_step );
				// Insert presteps
				Algorithm::poss_type
					pre_step_i = 0;
				// (10.-?) LineSearchFullStep
				algo->insert_assoc_step(
					step_num
					,GeneralIterationPack::PRE_STEP
					,++pre_step_i
					,"LineSearchFullStep"
					,line_search_full_step_step
					);
				// (10.-?) MeritFunc_PenaltyPramUpdate
				algo->insert_assoc_step(
					step_num
					,GeneralIterationPack::PRE_STEP
					,++pre_step_i
					,"MeritFunc_PenaltyParamUpdate"
					,merit_func_penalty_param_update_step
					);
				// ToDo: Add other LineSearch step objects
			}

		}
		else if ( mI > 0 || nb > 0 ) {
			//
			// Nonlinear inequality constrained NLP ( mI > 0 || num_bounded_x > 0 )
			//
			if(trase_out)
				*trase_out 
					<< "\nConfiguring an algorithm for nonlinear generally constrained "
					<< "NLP ( mI > 0 || num_bounded_x > 0 ) ...\n";

			THROW_EXCEPTION(
				mI > 0 || nb > 0, std::logic_error
				,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
				"General inequality constrained NLPS are not supported yet!" );
		}
		else {
			assert(0); // Error, this should not ever be called!
		}
	}
	
}

void rSQPAlgo_ConfigMamaJama::init_algo(rSQPAlgoInterface* _algo)
{
	using DynamicCastHelperPack::dyn_cast;
	namespace rcp = ReferenceCountingPack;

	THROW_EXCEPTION(
		_algo == NULL, std::invalid_argument
		,"rSQPAlgo_ConfigMamaJama::init_algo(_algo) : Error, "
		"_algo can not be NULL" );

	rSQPAlgo             &algo    = dyn_cast<rSQPAlgo>(*_algo);
	rSQPState	         &state   = algo.rsqp_state();
	NLP			         &nlp     = algo.nlp();
	NLPReduced           *nlp_red = dynamic_cast<NLPReduced*>(&nlp);
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
		const int num_opts = 15;
		enum EMamaJama {
			DIRECT_LINEAR_SOLVER
			,NULL_SPACE_MATRIX
			,RANGE_SPACE_MATRIX
			,MAX_BASIS_COND_CHANGE_FRAC
			,EXACT_REDUCED_HESSIAN
			,QUASI_NEWTON
			,MAX_DOF_QUASI_NEWTON_DENSE
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
			"direct_linear_solver"
			,"null_space_matrix"
			,"range_space_matrix"
			,"max_basis_cond_change_frac"
			,"exact_reduced_hessian"
			,"quasi_newton"
			,"max_dof_quasi_newton_dense"
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
				case DIRECT_LINEAR_SOLVER:
				{
					const std::string &linear_solver = ofsp::option_value(itr);
					if( linear_solver == "MA28" )
#ifdef SPARSE_SOLVER_PACK_USE_MA28
						ov->direct_linear_solver_type_ = LA_MA28;
#else
						throw std::logic_error(
							"rSQPAlgo_ConfigMamaJama::readin_options(...) : MA28 is not supported,"
							" must define SPARSE_SOLVER_PACK_USE_MA28!" );
#endif
					else if( linear_solver == "MA48" )
#ifdef SPARSE_SOLVER_PACK_USE_MA48
						ov->direct_linear_solver_type_ = LA_MA48;
#else
						throw std::logic_error(
							"rSQPAlgo_ConfigMamaJama::readin_options(...) : MA48 is not supported,"
							" must define SPARSE_SOLVER_PACK_USE_MA48!" );
#endif
					else if( linear_solver == "AUTO" )
						ov->direct_linear_solver_type_ = LA_AUTO;
					else
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"direct_linear_solver\" "
							"Only the options \'AUTO\' \'MA28\' and \'MA48\' are avalible." );
					break;
				}
				case NULL_SPACE_MATRIX:
				{
					const std::string &opt_val = ofsp::option_value(itr);
					if( opt_val == "EXPLICIT" )
						ov->null_space_matrix_type_ = NULL_SPACE_MATRIX_EXPLICIT;
					else if( opt_val == "IMPLICIT" )
						ov->null_space_matrix_type_ = NULL_SPACE_MATRIX_IMPLICIT;
					else if( opt_val == "AUTO" )
						ov->null_space_matrix_type_ = NULL_SPACE_MATRIX_AUTO;
					else
						throw std::invalid_argument( "crrSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"null_space_matrix\" "
							", Only the options for Z of EXPLICIT, IMPLICIT"
							", and AUTO are avalible."	);
					break;
				}
				case RANGE_SPACE_MATRIX:
				{
					const std::string &opt_val = ofsp::option_value(itr);
					if( opt_val == "COORDINATE" )
						ov->range_space_matrix_type_ = RANGE_SPACE_MATRIX_COORDINATE;
					else if( opt_val == "ORTHOGONAL" )
						ov->range_space_matrix_type_ = RANGE_SPACE_MATRIX_ORTHOGONAL;
					else if( opt_val == "AUTO" )
						ov->range_space_matrix_type_ = RANGE_SPACE_MATRIX_AUTO;
					else
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"range_space_matrix\" "
							", Only the options for Z of COORDINATE,"
							", ORTHOGONAL and AUTO are avalible."	);
					break;
				}
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
						throw std::invalid_argument( 
							"rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"quasi_newton\" "
							", Only options of BFGS, PBFGS"
							", LBFGS, LPBFGS and AUTO are avalible."
							);
					break;
				}
				case MAX_DOF_QUASI_NEWTON_DENSE:
					ov->max_dof_quasi_newton_dense_ = ::atoi( ofsp::option_value(itr).c_str() );
					break;
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
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"hessian_initialization\" "
							", Only options of IDENTITY, FINITE_DIFF_SCALE_IDENTITY,"
							" FINITE_DIFF_DIAGONAL, FINITE_DIFF_DIAGONAL_ABS and AUTO"
							" are available"  );
					break;
				}
				case QP_SOLVER:
				{
					const std::string &qp_solver = ofsp::option_value(itr);
					if( qp_solver == "AUTO" )
						ov->qp_solver_type_ = QP_AUTO;
					else if( qp_solver == "QPSOL" )
						ov->qp_solver_type_ = QP_QPSOL;
					else if( qp_solver == "QPOPT" )
#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
						ov->qp_solver_type_ = QP_QPOPT;
#else
						throw std::logic_error(
							"rSQPAlgo_ConfigMamaJama::readin_options(...) : QPOPT is not supported,"
							" must define CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT!" );
#endif
					else if( qp_solver == "QPKWIK" )
						ov->qp_solver_type_ = QP_QPKWIK;
					else if( qp_solver == "QPSCHUR" )
						ov->qp_solver_type_ = QP_QPSCHUR;
					else
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"qp_solver\" "
							"Only qp solvers QPOPT, QPSOL, QPKWIK, QPSCHUR and AUTO are avalible."	);
					break;
				}
				case REINIT_HESSIAN_ON_QP_FAIL:
					ov->reinit_hessian_on_qp_fail_ = StringToBool( "reinit_hessian_on_qp_fail", ofsp::option_value(itr).c_str() );
					break;
				case LINE_SEARCH_METHOD:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "NONE" )
						ov->line_search_method_ = LINE_SEARCH_NONE;
					else if( option == "DIRECT" )
						ov->line_search_method_ = LINE_SEARCH_DIRECT;
					else if( option == "2ND_ORDER_CORRECT" )
						ov->line_search_method_ = LINE_SEARCH_2ND_ORDER_CORRECT;
					else if( option == "WATCHDOG" )
						ov->line_search_method_ = LINE_SEARCH_WATCHDOG;
					else if( option == "AUTO" )
						ov->line_search_method_ = LINE_SEARCH_AUTO;
					else
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
							"Error, incorrect value for \"line_search_method\".\n"
							"Only the options NONE, DIRECT, 2ND_ORDER_CORRECT, WATCHDOG "
							"and AUTO are avalible." );
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
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
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
						throw std::invalid_argument( "rSQPAlgo_ConfigMamaJama::readin_options(...) : "
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

	if( cov->direct_linear_solver_type_ == LA_AUTO && uov.direct_linear_solver_type_ == LA_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\ndirect_linear_solver == AUTO: setting direct_linear_solver = MA28\n";
		cov->direct_linear_solver_type_ = LA_MA28;
	}
	else if( cov->direct_linear_solver_type_ == LA_AUTO ) {
		cov->direct_linear_solver_type_ = uov.direct_linear_solver_type_;
	}
	if( cov->null_space_matrix_type_ == NULL_SPACE_MATRIX_AUTO && uov.null_space_matrix_type_ == NULL_SPACE_MATRIX_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nnull_space_matrix_type == AUTO: Let the algorithm deside as it goes along\n";
	}
	else if(cov->null_space_matrix_type_ == NULL_SPACE_MATRIX_AUTO) {
		cov->null_space_matrix_type_ = uov.null_space_matrix_type_;
	}
	if( cov->range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO && uov.range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nrange_space_matrix_type == AUTO: Let the algorithm deside as it goes along\n";
	}
	else if(cov->range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO) {
		cov->range_space_matrix_type_ = uov.range_space_matrix_type_;
	}
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
	if( cov->max_dof_quasi_newton_dense_ < 0 && uov.max_dof_quasi_newton_dense_ < 0 ) {
		if(trase_out)
			*trase_out
				<< "\nmax_dof_quasi_newton_dense < 0 : setting max_dof_quasi_newton_dense = 500\n";
		cov->max_dof_quasi_newton_dense_ = 500;
	}
	else if(cov->max_dof_quasi_newton_dense_ < 0) {
		cov->max_dof_quasi_newton_dense_ = uov.max_dof_quasi_newton_dense_;
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
				<< "\nhessian_initialization == AUTO: setting hessian_initialization = FINITE_DIFF_DIAGONAL_ABS\n";
		cov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS;
	}
	else if(cov->hessian_initialization_ == INIT_HESS_AUTO) {
		cov->hessian_initialization_ = uov.hessian_initialization_;
	}
	if( cov->qp_solver_type_ == QP_AUTO && uov.qp_solver_type_ == QP_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nqp_solver_type == AUTO: setting qp_solver_type = QPSCHUR\n";
		cov->qp_solver_type_ = QP_QPSCHUR;
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
