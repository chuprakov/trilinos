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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)

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

//#include "ReducedSpaceSQPPack/include/std/DecompositionSystemVarReductStd.h"
//#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateDirect.h"
//#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateAdjoint.h"

//#include "SparseSolverPack/test/BasisSystemTesterSetOptions.h"

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

//#include "ReducedSpaceSQPPack/include/std/DecompositionSystemVarReductStd.h"
//#include "ReducedSpaceSQPPack/include/std/EvalNewPointStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachCoordinate_Step.h"
//#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachOrthogonal_Step.h"
//#include "ReducedSpaceSQPPack/include/std/ReducedGradientStd_Step.h"
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
//#include "ReducedSpaceSQPPack/include/std/DepDirecStd_Step.h"
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

//#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateDirect.h"
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
// should agree with what are in the rSQPpp.opt.rsqp_mama_jama_solve file.
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
	const basis_sys_ptr_t&  basis_sys
	)
	:basis_sys_(basis_sys)
{}

rSQPAlgo_ConfigMamaJama::~rSQPAlgo_ConfigMamaJama() {
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
			<< "*********************************************\n"
			<< "*** rSQPAlgo_ConfigMamaJama configuration ***\n"
			<< "*********************************************\n";
	}

	// ////////////////////////////////////////////////////////////
	// A. ???

	// /////////////////////////////////////////////////////////////////////////
	// B. Create an algo object, give to algo_cntr, then give algo_cntr to algo

	if(trase_out)
		*trase_out << "\nCreating the algo object ...\n";

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
		*trase_out << "\nSetting the NLP and track objects ...\n";

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
		tailored_approach = false;
	}
	else {
		if( nlp_fod ) {
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

	// //////////////////////////////////////////////////////
	// C.1.  Sort out the options

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

	if( uov_.range_space_matrix_type_ != RANGE_SPACE_MATRIX_COORDINATE ) {
		if(trase_out)
			*trase_out <<
				"\nrange_space_matrix != COORDINATE:\n"
				"Sorry, the orthogonal decomposition is not updated yet!\n"
				"setting range_space_matrix = COORDINATE ...\n";
		cov_.range_space_matrix_type_ = RANGE_SPACE_MATRIX_COORDINATE;
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
	// C.1. Create and set the state object

	if(trase_out)
		*trase_out << "\nCreating and setting the state object ...\n";

	{
		//
		// Create the state object with the vector spaces
		//

		typedef ref_count_ptr<rSQPState>   state_ptr_t;
		state_ptr_t state = NULL;
		if( tailored_approach ) {
			state = rcp::rcp(
				new rSQPState(
					nlp.space_x()
					,nlp.space_c()
					,nlp.space_h()
					,( nlp_fod->var_dep().size() 
					   ? nlp.space_x()->sub_space(nlp_fod->var_dep())->clone()
					   : NULL )
					,( nlp_fod->var_indep().size()
					   ? nlp.space_x()->sub_space(nlp_fod->var_indep())->clone()
					   : NULL )
					)
				);
		}
		else {
			assert(0); // ToDo: Implement by adding space_py() and space_pz() to BasisSystem!
		}
		
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

		// ToDo: Add matrix iq object for R
		// ToDo: Add matrix iq object for Uz
		// ToDo: Add matrix iq object for Uy
		// ToDo: Add matrix iq object for Vz
		// ToDo: Add matrix iq object for Vy
		// ToDo: Add other needed projected matrix iq objects

		// Add reduced Hessian

		if( !cov_.exact_reduced_hessian_ ) {
			ref_count_ptr<afp::AbstractFactory<MatrixSymWithOp> >
				abstract_factory_rHL = NULL;
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
				merit_func_factory = NULL;
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
	// C.2. Create and set the decomposition system object

	if(!tailored_approach) {
		assert(0); // ToDo: Create the proper decomp_sys object
	}

	// /////////////////////////////////////////////////////
	// C.3  Create and set the step objects

	if(trase_out)
		*trase_out << "\nCreating and setting the step objects ...\n";

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
		bounds_tester_ptr_t
			bounds_tester = NULL;
		if(nb) { // has variable bounds?
			const value_type var_bounds_warning_tol = 1e-10;
			bounds_tester = rcp::rcp(
				new VariableBoundsTester(
					var_bounds_warning_tol              // default warning tolerance
					, algo_cntr->max_var_bounds_viol()	// default warning tolerance
					) );
			if(options_.get()) {
				ConstrainedOptimizationPack::VariableBoundsTesterSetOptions
					options_setter( bounds_tester.get() );
				options_setter.set_options(*options_);
				}
		}

		// EvalNewPoint_Step
		algo_step_ptr_t    eval_new_point_step = NULL;
		{
			// Create the step object
			if( tailored_approach ) {
				// create and setup the derivative tester
				typedef NLPInterfacePack::NLPFirstOrderDirectTester     NLPFirstOrderDirectTester;
				typedef rcp::ref_count_ptr<NLPFirstOrderDirectTester>   deriv_tester_ptr_t;
				deriv_tester_ptr_t
					deriv_tester = rcp::rcp(
						new NLPFirstOrderDirectTester(
							NLPFirstOrderDirectTester::FD_DIRECTIONAL    // Gf testing
							,NLPFirstOrderDirectTester::FD_DIRECTIONAL   // -Inv(C)*N testing
							) );
				if(options_.get()) {
					NLPInterfacePack::NLPFirstOrderDirectTesterSetOptions
						options_setter(deriv_tester.get());
					options_setter.set_options(*options_);
				}
				// create the step
				typedef rcp::ref_count_ptr<EvalNewPointTailoredApproach_Step>  _eval_new_point_step_ptr_t;
				_eval_new_point_step_ptr_t
					_eval_new_point_step = NULL;
				switch( cov_.range_space_matrix_type_ ) {
					case RANGE_SPACE_MATRIX_COORDINATE:
						_eval_new_point_step
							= rcp::rcp(new EvalNewPointTailoredApproachCoordinate_Step(deriv_tester,bounds_tester));
						break;
					case RANGE_SPACE_MATRIX_ORTHOGONAL:
//						eval_new_point_step
//							= new EvalNewPointTailoredApproachOrthogonal_Step(deriv_tester,bounds_tester);
						THROW_EXCEPTION(
							true, std::logic_error
							,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
							"the othogonal decomposition is not updated for NLPFirstOrderDirect yet!" );
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
				THROW_EXCEPTION(
					true, std::logic_error
					,"rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : Error, "
					"NLPFirstOrderInfo is not supported yet!" );
			}
		}

		// RangeSpace_Step
		algo_step_ptr_t    range_space_step_step = NULL;
		if( !tailored_approach ) {
			assert(0); // ToDo: Implment!
		}

		// ReducedGradient_Step
		algo_step_ptr_t    reduced_gradient_step = NULL;
		if( !tailored_approach ) {
			assert(0); // ToDo: Implment!
		}

		// CheckSkipBFGSUpdate
		algo_step_ptr_t    check_skip_bfgs_update_step = NULL;
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
		algo_step_ptr_t    reduced_hessian_step = NULL;
		{
			// Get the strategy object that will perform the actual secant update.
			ref_count_ptr<ReducedHessianSecantUpdate_Strategy>
				secant_update_strategy = NULL;
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
				,new IterQuantityAccessContiguous<QuasiNewtonStats>(
					1, quasi_newton_stats_name )
				);
		}

		// NullSpace_Step
		algo_step_ptr_t    null_space_step_step = NULL;
		{
			null_space_step_step = rcp::rcp(new NullSpaceStepWithoutBounds_Step());
		}

		// CalcDFromYPYZPZ_Step
		algo_step_ptr_t    calc_d_from_Ypy_Zpy_step = NULL;
		{
			calc_d_from_Ypy_Zpy_step = rcp::rcp(new CalcDFromYPYZPZ_Step());
		}

		// CalcReducedGradLagrangianStd_AddedStep
		algo_step_ptr_t    calc_reduced_grad_lagr_step = NULL;
		{
			calc_reduced_grad_lagr_step = rcp::rcp(
				new CalcReducedGradLagrangianStd_AddedStep() );
		}

		// CheckConvergence_Step
		algo_step_ptr_t    check_convergence_step = NULL;
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
		algo_step_ptr_t    merit_func_penalty_param_update_step = NULL;
		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<MeritFunc_PenaltyParamUpdate_AddedStep>
				param_update_step = NULL;
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
		algo_step_ptr_t    line_search_full_step_step = NULL;
		{
			line_search_full_step_step = rcp::rcp(new LineSearchFullStep_Step(bounds_tester));
		}

		// LineSearch_Step
		algo_step_ptr_t    line_search_step = NULL;
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
				THROW_EXCEPTION(
					m == 0 && mI == 0 && nb == 0, std::logic_error
					,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
					"Unconstrained NLPs are not supported yet!" );
			}
			else {
				//
				// Simple bound constrained NLP (m == 0, mI == 0, num_bounded_x > 0)
				//
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

			int step_num = 0;
	
			// (1) EvalNewPoint
			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );

			// (2) CalcReducedGradLagrangian
			algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

			// (3) CalcLagrangeMultDecomposed
			// Compute these here so that in case we converge we can report them
			if( !tailored_approach ) {
				assert(0); // ToDo: Insert this step
			}

			// (4) CheckConvergence
			algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );

			// (5) ReducedHessian
			algo->insert_step( ++step_num, ReducedHessian_name, reduced_hessian_step );

			// (5.-1) CheckSkipBFGSUpdate
			algo->insert_assoc_step(
				step_num
				,GeneralIterationPack::PRE_STEP
				,1
				,CheckSkipBFGSUpdate_name
				,check_skip_bfgs_update_step
			  );

			// (6) NullSpaceStep
			algo->insert_step( ++step_num, NullSpaceStep_name, null_space_step_step );

			// (7) CalcDFromYPYZPZ
			algo->insert_step( ++step_num, CalcDFromYPYZPZ_name, calc_d_from_Ypy_Zpy_step );
			// ToDo: Add ths step!
			
			// (8) LineSearch
			if( cov_.line_search_method_ == LINE_SEARCH_NONE ) {
				algo->insert_step( ++step_num, LineSearch_name, line_search_full_step_step );
			}
			else {
				// (8) Main line search step
				algo->insert_step( ++step_num, LineSearch_name, line_search_step );
				// Insert presteps
				Algorithm::poss_type
					pre_step_i = 0;
				// (8.-?) LineSearchFullStep
				algo->insert_assoc_step(
					step_num
					,GeneralIterationPack::PRE_STEP
					,++pre_step_i
					,"LineSearchFullStep"
					,line_search_full_step_step
					);
				// (8.-?) MeritFunc_PenaltyPramUpdate
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
			THROW_EXCEPTION(
				mI > 0 || nb > 0, std::logic_error
				,"rSQPAlgo_ConfigMamaJama::config_alg_cntr(...) : Error, "
				"General inequality constrained NLPS are not supported yet!" );
		}
		else {
			assert(0); // Error, this should not ever be called!
		}

		// Set the standard set of step objects.
		// This default algorithm is for NLPs with no variable bounds.
/*

		// (1) EvalNewPoint
		if( tailored_approach ) {
			typedef rcp::ref_count_ptr<NLPrSQPTailoredApproachTester>
				nonconst_deriv_tester_ptr_t;
			typedef EvalNewPointTailoredApproach_Step::deriv_tester_ptr_t
				deriv_tester_ptr_t;

			nonconst_deriv_tester_ptr_t
				deriv_tester = new NLPrSQPTailoredApproachTester(
										NLPrSQPTailoredApproachTester::FD_DIRECTIONAL		// Gf
										,NLPrSQPTailoredApproachTester::FD_COMPUTE_ALL		// -Inv(C)*N
										);

			{
				NLPrSQPTailoredApproachTesterSetOptions options_setter(deriv_tester.get());
				if(options_) options_setter.set_options(*options_);
			}		

			EvalNewPointTailoredApproach_Step
				*eval_new_point_step = NULL;
			switch( cov_.range_space_matrix_type_ ) {
				case RANGE_SPACE_MATRIX_COORDINATE:
					eval_new_point_step
						= new EvalNewPointTailoredApproachCoordinate_Step(deriv_tester,bounds_tester);
					break;
				case RANGE_SPACE_MATRIX_ORTHOGONAL:
					eval_new_point_step
						= new EvalNewPointTailoredApproachOrthogonal_Step(deriv_tester,bounds_tester);
					break;
				default:
					assert(0);	// only a local error
			}
			{
				EvalNewPointTailoredApproach_StepSetOptions options_setter(eval_new_point_step);
				if(options_) options_setter.set_options(*options_);
			}		

			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );
		}
		else {
			typedef NLPInterfacePack::TestingPack::NLPFirstDerivativesTester
				NLPFirstDerivativesTester;
			typedef NLPInterfacePack::TestingPack::NLPFirstDerivativesTesterSetOptions
				NLPFirstDerivativesTesterSetOptions;
			typedef rcp::ref_count_ptr<NLPFirstDerivativesTester>
				nonconst_deriv_tester_ptr_t;
			typedef EvalNewPointStd_Step::deriv_tester_ptr_t
				deriv_tester_ptr_t;
			
			nonconst_deriv_tester_ptr_t
				deriv_tester = new NLPFirstDerivativesTester();

			{
				NLPFirstDerivativesTesterSetOptions options_setter(deriv_tester.get());
				if(options_) options_setter.set_options(*options_);
			}		

			EvalNewPointStd_Step
				*eval_new_point_step = new EvalNewPointStd_Step(deriv_tester,bounds_tester);

			{
				EvalNewPointStd_StepSetOptions options_setter(eval_new_point_step);
				if(options_) options_setter.set_options(*options_);
			}		

			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );
		}

		// (2) DepDirec
		if( !tailored_approach ) {
			algo->insert_step( ++step_num, DepDirec_name, new  DepDirecStd_Step );
		}

		// (3) ReducedGradient
		algo->insert_step( ++step_num, ReducedGradient_name, new ReducedGradientStd_Step );

		// (4) Calculate Reduced Gradient of the Lagrangian
		algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, new CalcReducedGradLagrangianStd_AddedStep );

		// (5)	Calculate the Lagrange multipliers for the independent constraints.
		// 		These are computed here just in case the algorithm converges and we need to
		// 		report these multipliers to the NLP.
		if( !tailored_approach ) {
			algo->insert_step( ++step_num, CalcLambdaIndep_name, new  CalcLambdaIndepStd_AddedStep );
		}

		// (6) Check for convergence
		{
			CheckConvergenceStd_AddedStep
				*check_convergence_step = new CheckConvergenceStd_AddedStep;

			CheckConvergenceStd_AddedStepSetOptions
				opt_setter( check_convergence_step );
			if(options_) opt_setter.set_options( *options_ );

			algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
		}

		// (7) ReducedHessian
		{
			// Get the strategy object that will perform the actual secant update.
			ref_count_ptr<ReducedHessianSecantUpdate_Strategy>
				secant_update_strategy = NULL;
			switch( cov_.quasi_newton_ )
			{
			    case QN_BFGS:
			    case QN_PBFGS:
			    case QN_LBFGS:
			    case QN_LPBFGS:
				{
					typedef ref_count_ptr<BFGSUpdate_Strategy> bfgs_strategy_ptr_t;
					bfgs_strategy_ptr_t
						bfgs_strategy = new BFGSUpdate_Strategy;
					BFGSUpdate_StrategySetOptions
						opt_setter( bfgs_strategy.get() );
					if(options_) opt_setter.set_options( *options_ );

					switch( cov_.quasi_newton_ ) {
					    case QN_BFGS:
					    case QN_LBFGS:
						{
							secant_update_strategy = new ReducedHessianSecantUpdateBFGSFull_Strategy(bfgs_strategy);
							break;
						}
					    case QN_PBFGS:
					    case QN_LPBFGS:
						{
							typedef ref_count_ptr<ReducedHessianSecantUpdateBFGSProjected_Strategy>
								pbfgs_strategy_ptr_t;
							pbfgs_strategy_ptr_t
								pbfgs_strategy = new ReducedHessianSecantUpdateBFGSProjected_Strategy(bfgs_strategy);
							ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions
								opt_setter( pbfgs_strategy.get() );
							if(options_) opt_setter.set_options( *options_ );
							if( cov_.quasi_newton_ == QN_PBFGS ) {
								secant_update_strategy
									= rcp::rcp_static_cast<ReducedHessianSecantUpdate_Strategy>(pbfgs_strategy);
							}
							else {
								typedef ref_count_ptr<ReducedHessianSecantUpdateLPBFGS_Strategy>
									lpbfgs_strategy_ptr_t;
								lpbfgs_strategy_ptr_t
									lpbfgs_strategy
									= new ReducedHessianSecantUpdateLPBFGS_Strategy(
										pbfgs_strategy 
										);
								ReducedHessianSecantUpdateLPBFGS_StrategySetOptions
									opt_setter( lpbfgs_strategy.get() );
								if(options_) opt_setter.set_options( *options_ );
								secant_update_strategy
									= rcp::rcp_static_cast<ReducedHessianSecantUpdate_Strategy>(lpbfgs_strategy);
							}
							break;
						}
					}
					break;
				}
			    default:
					assert(0);
			}
			// Finally build the step object
			ReducedHessianSecantUpdateStd_Step
				*step = new ReducedHessianSecantUpdateStd_Step( secant_update_strategy );
			algo->insert_step( ++step_num, ReducedHessian_name, step );
		}

		// (7.-1) CheckSkipBFGSUpdate
		{
			CheckSkipBFGSUpdateStd_Step
				*step = new CheckSkipBFGSUpdateStd_Step;
			CheckSkipBFGSUpdateStd_StepSetOptions
				opt_setter( step );
			if(options_) opt_setter.set_options( *options_ );
			
			algo->insert_assoc_step(
				  step_num
				, GeneralIterationPack::PRE_STEP
				, 1
				, CheckSkipBFGSUpdate_name
				, step
			  );
		}

		// (8) IndepDirec
		algo->insert_step( ++step_num, IndepDirec_name, new  IndepDirecWithoutBounds_Step );

		// (9) CalcDFromYPYZPZ
		algo->insert_step( ++step_num, "CalcDFromYPYZPZ", new CalcDFromYPYZPZ_Step );

		// (10) LineSearch
		if( cov_.line_search_method_ == LINE_SEARCH_NONE ) {

			algo->insert_step(
				  ++step_num
				, LineSearch_name
				, new LineSearchFullStep_Step(bounds_tester)
			  );
			
		}
		else {

			using ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy;
			using ConstrainedOptimizationPack::MeritFuncNLPL1;
			using ConstrainedOptimizationPack::MeritFuncNLPModL1;

			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1:
					merit_func = merit_func_ptr_t(new MeritFuncNLPL1);
					break;
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
					merit_func = merit_func_ptr_t(new MeritFuncNLPModL1);
					break;
				default:
					assert(0);	// local programming error
			}

			DirectLineSearchArmQuad_Strategy
				*direct_line_search = new  DirectLineSearchArmQuad_Strategy();

			ConstrainedOptimizationPack::DirectLineSearchArmQuad_StrategySetOptions
				ls_options_setter( direct_line_search, "DirectLineSearchArmQuadSQPStep" );
			if(options_) ls_options_setter.set_options( *options_ );
			LineSearch_Step
				*line_search_step = 0;

			switch( cov_.line_search_method_ ) {
				case LINE_SEARCH_DIRECT: {
					line_search_step = new LineSearchDirect_Step(
						  direct_line_search
						, rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
					  );
					break;
				}
				case LINE_SEARCH_2ND_ORDER_CORRECT: {
					typedef LineSearch2ndOrderCorrect_Step ls_t;
					// Set the default print level for the newton iterations
					ls_t::ENewtonOutputLevel
						newton_olevel;
					switch( algo_cntr->journal_output_level() ) {
						case PRINT_NOTHING:
							newton_olevel = ls_t::PRINT_NEWTON_NOTHING;
							break;
						case PRINT_BASIC_ALGORITHM_INFO:
							newton_olevel = ls_t::PRINT_NEWTON_NOTHING;
							break;
						case PRINT_ALGORITHM_STEPS:
						case PRINT_ACTIVE_SET:
							newton_olevel = ls_t::PRINT_NEWTON_SUMMARY_INFO;
							break;
						case PRINT_VECTORS:
						case PRINT_ITERATION_QUANTITIES:
							newton_olevel = ls_t::PRINT_NEWTON_VECTORS;
							break;
					}
					// Create the QuasiRangeSpaceStep_Strategy object
					QuasiRangeSpaceStep_Strategy
						*quasi_range_space_step = NULL;
					if(tailored_approach) {
						// ToDo: Set the EvalNewPointTailoredApproach_Step object!
						quasi_range_space_step = new QuasiRangeSpaceStepTailoredApproach_Strategy();
					}
					else {
						quasi_range_space_step = new QuasiRangeSpaceStepStd_Strategy();
					}
					// Create the FeasibilityStepReducedStd_Strategy object and
					// set its options
					{
						feasibility_step_strategy
							= new FeasibilityStepReducedStd_Strategy(
								quasi_range_space_step         // Given ownership to delete!
								,NULL  // QP solver (must be set later)
								,NULL  // QP solver tester (must be set later)
								);
						// Set the options
						FeasibilityStepReducedStd_StrategySetOptions
							opt_setter(feasibility_step_strategy.get());
						if(options_) opt_setter.set_options( *options_ );
					}
					// Create the line search object for the newton iterations
					// and set its options from the option from stream object.
					DirectLineSearchArmQuad_Strategy
						*direct_ls_newton = new  DirectLineSearchArmQuad_Strategy(
								5	// max_iter
							);
					ConstrainedOptimizationPack::DirectLineSearchArmQuad_StrategySetOptions
						direct_ls_opt_setter( direct_ls_newton
							, "DirectLineSearchArmQuad2ndOrderCorrectNewton" );
					if(options_) direct_ls_opt_setter.set_options( *options_ );
					// Create the step object and set its options from the options object.
					LineSearch2ndOrderCorrect_Step
						*_line_search_step = new LineSearch2ndOrderCorrect_Step(
							direct_line_search
							,rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
							,rcp::rcp_implicit_cast<FeasibilityStep_Strategy>(
								feasibility_step_strategy)
							,direct_ls_newton
							,direct_line_search->eta()
							,newton_olevel
						);
					LineSearch2ndOrderCorrect_StepSetOptions
						ls_opt_setter( _line_search_step );
					if(options_) ls_opt_setter.set_options( *options_ );
					//
					line_search_step = _line_search_step;
					break;
				}
				case LINE_SEARCH_WATCHDOG: {
					// Create the step object and set its options from the options object.
					LineSearchWatchDog_Step
						*_line_search_step = new LineSearchWatchDog_Step(
							  direct_line_search
							, rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
							, direct_line_search->eta()
							);
					LineSearchWatchDog_StepSetOptions
						ls_opt_setter( _line_search_step );
					if(options_) ls_opt_setter.set_options( *options_ );
					//
					line_search_step = _line_search_step;
					break;
				}
			}

			algo->insert_step(
				++step_num
				, LineSearch_name
				, new LineSearchFailureNewBasisSelection_Step( 
					line_search_step
					,new NewBasisSelectionStd_Strategy(decomp_sys.get()) // Just aggregation, okay!
					)
				);
			
			Algorithm::poss_type
				pre_step_i = 0;

			// (10.-?) Update the penalty parameter for the Merit function.
			MeritFunc_PenaltyParamUpdate_AddedStep
				*param_update_step = 0;

			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1: {
					switch(cov_.l1_penalty_param_update_) {
						case L1_PENALTY_PARAM_WITH_MULT:
							param_update_step
								= new  MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
											rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
							break;
						case L1_PENALTY_PARAM_MULT_FREE:
							param_update_step
								= new  MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
											rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
							break;
						default:
							assert(0);
					}
					break;
				}
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
					param_update_step = new  MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(
											rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
					break;
				default:
					assert(0);	// local programming error
			}
			// Set the options from the stream
			MeritFunc_PenaltyParamUpdate_AddedStepSetOptions
				ppu_options_setter( param_update_step );
			if(options_) ppu_options_setter.set_options( *options_ );

			algo->insert_assoc_step(	  step_num
										, GeneralIterationPack::PRE_STEP
										, ++pre_step_i
										, "MeritFunc_PenaltyParamUpdate"
										, param_update_step // give control over memory!
									);
			
			// (10.-?)	Compute the full step before the linesearch
			algo->insert_assoc_step(	  step_num
										, GeneralIterationPack::PRE_STEP
										, ++pre_step_i
										, "LineSearchFullStep"
										, new  LineSearchFullStep_Step(bounds_tester)	);

			// (10.-?) Increase the penalty parameters to get a larger step.
			if( cov_.merit_function_type_ == MERIT_FUNC_MOD_L1_INCR ) {

				MeritFunc_ModifiedL1LargerSteps_AddedStep
					*_added_step = new MeritFunc_ModifiedL1LargerSteps_AddedStep(
						   rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
						 , direct_line_search->eta() );

				MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions
					options_setter( _added_step );
				if(options_) options_setter.set_options( *options_ );

				algo->insert_assoc_step(	  step_num
											, GeneralIterationPack::PRE_STEP
											, ++pre_step_i
											, "MeritFunc_ModifiedL1LargerSteps"
											, _added_step // give control over memory!
										);
			}
		}
	
*/

/*
		// Reconfigure the steps if the NLP has bounds

		if( algo_cntr->nlp().has_bounds() ) {

			if(trase_out)
				*trase_out << "\nReconfiguring steps for NLP with bounds ...\n";
			
			// If the NLP has bounds on the variables then we must replace
			// IndepDirec and move the convergence check to the end (along
			// with the calculation of the reduced gradient of the lagrangian).
			
			// Setup IndepDirec step
			
			// Add iteration quantity for d_bounds
			algo->state().set_iter_quant( d_bounds_name
				, new IterQuantityAccessContinuous<SparseBounds>( 1, d_bounds_name ) );

			// Create the QP solver

			// RAB 8/28/00: In the future, all the QP solvers will also use this interface.

			switch( cov_.qp_solver_type_ ) {
				case QP_QPSCHUR: {
					using ConstrainedOptimizationPack::QPSolverRelaxedQPSchur;
					using ConstrainedOptimizationPack::QPSolverRelaxedQPSchurSetOptions;
					using ConstrainedOptimizationPack::QPSchurInitKKTSystemHessianFull;
					using ConstrainedOptimizationPack::QPSchurInitKKTSystemHessianSuperBasic;
					ConstrainedOptimizationPack::QPSolverRelaxedQPSchur::InitKKTSystem
						*init_kkt_sys = NULL;
					switch( cov_.quasi_newton_ ) {
						case QN_BFGS:
						case QN_LBFGS:
							init_kkt_sys = new QPSchurInitKKTSystemHessianFull();
							break;
						case QN_PBFGS:
						case QN_LPBFGS:
							init_kkt_sys = new QPSchurInitKKTSystemHessianSuperBasic();
							break;
						default:
							assert(0);
					}
					QPSolverRelaxedQPSchur
						*_qp_solver = new QPSolverRelaxedQPSchur(init_kkt_sys);
					QPSolverRelaxedQPSchurSetOptions
						qp_options_setter(_qp_solver);
					qp_options_setter.set_options( *options_ );
					qp_solver = _qp_solver; // give ownership to delete!
					break;
				}
				case QP_QPKWIK: {
					using ConstrainedOptimizationPack::QPSolverRelaxedQPKWIK;
					QPSolverRelaxedQPKWIK
						*_qp_solver = new QPSolverRelaxedQPKWIK();
					qp_solver = _qp_solver; // give ownership to delete!
					break;
				}
				case QP_QPOPT: {
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
					break;
				}
				default:
					assert(0);
			}

			// Create the QP solver tester object and set its options
			qp_tester = new QPSolverRelaxedTester();
			ConstrainedOptimizationPack::QPSolverRelaxedTesterSetOptions
				qp_tester_options_setter( qp_tester.get() );
			qp_tester_options_setter.set_options( *options_ );

			// Create the step object

			IndepDirecWithBoundsStd_Step
				*indep_direct_step = new  IndepDirecWithBoundsStd_Step( qp_solver, qp_tester );
			IndepDirecWithBoundsStd_StepSetOptions
				indep_direct_step_options_setter( indep_direct_step );
			indep_direct_step_options_setter.set_options( *options_ );

			// Set the step objects

			if( trase_out ) {
				*trase_out
					<< "\n\nreinit_hessian_on_qp_fail = " << (cov_.reinit_hessian_on_qp_fail_ ? "true" : "false");
				if(cov_.reinit_hessian_on_qp_fail_)
					*trase_out
						<< "\nThe algorithm will reinitalize the reduced hessian if the QP subproblem fails"
						<< "\nand will not throw a QPFailure exception the first time\n";
				else
					*trase_out
						<< "\nThe algorithm will not reinitalize the reduced hessian if the QP subproblem fails"
						<< "\nIt will just be a QPFailure exception thrown\n";
			}
			rSQPAlgo_Step
				*_indep_direc_step = NULL;
			if(cov_.reinit_hessian_on_qp_fail_)
				_indep_direc_step = new QPFailureReinitReducedHessian_Step(indep_direct_step); // Will clean up memory!
			else
				_indep_direc_step = indep_direct_step;
			Algorithm::poss_type poss;
			poss = algo->get_step_poss(IndepDirec_name);
			algo->remove_step( poss );	// remove any pre or post steps also
			algo->insert_step(
				poss
				,IndepDirec_name
				,_indep_direc_step // Will clean up memory!
				);
			algo->insert_assoc_step(
				 poss
				, GeneralIterationPack::PRE_STEP
				, 1
				, "SetDBounds"
				, new  SetDBoundsStd_AddedStep
				);

			Algorithm::step_ptr_t calc_rgrad_lagr, calc_lambda, check_conv;

			// Remove and save CalcReducedGradLagr..., CalcLambdaIndep... and CheckConv...

			poss = algo->get_step_poss(CalcReducedGradLagrangian_name);
			calc_rgrad_lagr	= algo->get_step(poss);
			algo->remove_step(poss);

			if(!tailored_approach) {
				poss			= algo->get_step_poss(CalcLambdaIndep_name);
				calc_lambda		= algo->get_step(poss);
				algo->remove_step(poss);
			}

			poss			= algo->get_step_poss(CheckConvergence_name);
			check_conv		= algo->get_step(poss);
			algo->remove_step(poss);

			// Add them before LineSearch
	
			poss		= algo->get_step_poss(LineSearch_name);
			algo->insert_step( poss++, CalcReducedGradLagrangian_name, calc_rgrad_lagr );
			if(!tailored_approach) {
				algo->insert_step( poss++, CalcLambdaIndep_name, calc_lambda );
			}
			algo->insert_step( poss++, CheckConvergence_name, check_conv );

		}

*/

	}

/*

	// 8/30/99: Add the iteration quantities for the QPSolverStats and 
	// ActSetStats and the added step that will calculate the
	// active set changes.

	{
		// Add active set and QP statistics to state object
		algo->state().set_iter_quant( act_set_stats_name
			, new IterQuantityAccessContinuous<ActSetStats>( 1, act_set_stats_name ) );
		algo->state().set_iter_quant( qp_solver_stats_name
			, new IterQuantityAccessContinuous<QPSolverStats>( 1, qp_solver_stats_name ) );

		if( algo_cntr->nlp().has_bounds() ) {

			// Insert active set computational step into algorithm
			Algorithm::poss_type
				poss = algo->get_step_poss( IndepDirec_name );
			if( poss == Algorithm::DOES_NOT_EXIST ) {
				// If we are not using IndepDirec to compute the active set then
				// we must be using SearchDirec
				assert( poss = algo->get_step_poss( SearchDirec_name ) );		
			}
			Algorithm::poss_type
				nas = algo->num_assoc_steps( poss, GeneralIterationPack::POST_STEP );
			algo->insert_assoc_step( poss, GeneralIterationPack::POST_STEP, nas+1
				, "ActiveSetStatistics", new ActSetStats_AddedStep );
			
			// Output the number of fixed depenent and indepent variables.
			algo->insert_assoc_step( poss, GeneralIterationPack::POST_STEP, nas+2
				, "CountNumFixedDepIndep", new NumFixedDepIndep_AddedStep );
		}

	}

	// 10/15/99: Add basis checking step
	if( !tailored_approach && cov_.max_basis_cond_change_frac_ < INF_BASIS_COND_CHANGE_FRAC ) {
		// 3/7/01: Added another test also
		CheckBasisFromCPy_Step
			*basis_check_step1 = new CheckBasisFromCPy_Step(
				new NewBasisSelectionStd_Strategy(decomp_sys.get()) );
		CheckBasisFromPy_Step
			*basis_check_step2 = new CheckBasisFromPy_Step(
				new NewBasisSelectionStd_Strategy(decomp_sys.get()) );
		if( cov_.max_basis_cond_change_frac_ > 0.0 ) {
			basis_check_step1->max_basis_cond_change_frac( cov_.max_basis_cond_change_frac_ );
			basis_check_step2->max_basis_cond_change_frac( cov_.max_basis_cond_change_frac_ );
		}
		Algorithm::poss_type poss;
		assert(poss = algo->get_step_poss( DepDirec_name ) );
		algo->insert_step(
			  ++poss
			, "CheckBasis1"
			, basis_check_step1
		  );
		algo->insert_step(
			  ++poss
			, "CheckBasis2"
			, basis_check_step2
		  );
	}

	// 10/19/99: Add quasi-newton stats
	algo->state().set_iter_quant( quasi_newton_stats_name
		, new IterQuantityAccessContinuous<QuasiNewtonStats>( 1, quasi_newton_stats_name ) );

	// 11/12/99: Add the option of using full steps after a specified number of iterations
	// if this option has been set and if we are using a linesearch method.
	if( cov_.full_steps_after_k_ > 0 && cov_.line_search_method_ != LINE_SEARCH_NONE ) {
		Algorithm::poss_type poss;
		assert(poss = algo->get_step_poss( LineSearch_name ) );
		Algorithm::step_ptr_t
			existing_step_ptr = algo->get_step( poss );
		// Check that the existing step is indeed a linesearch step and that
		// we can safely use static cast cast up in types.  If multiple
		// inheritance was used then this test will also fail.
		assert( dynamic_cast<LineSearch_Step*>(existing_step_ptr.get())
				==  existing_step_ptr.get() );
		rcp::ref_count_ptr<LineSearch_Step>
			existing_ls_step_ptr = rcp::rcp_static_cast<LineSearch_Step>(existing_step_ptr);
		algo->replace_step(
			  poss
			, new LineSearchFullStepAfterKIter_Step(
				  existing_ls_step_ptr
				, cov_.full_steps_after_k_
			  )
		  );

	}

	// 12/3/99: Adding finite difference initializaiton of the reduced hessian
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
		InitFinDiffReducedHessian_Step
			*init_red_hess_step = new InitFinDiffReducedHessian_Step(init_hess);
		InitFinDiffReducedHessian_StepSetOptions
			opt_setter( init_red_hess_step );
		if(options_) opt_setter.set_options( *options_ );
		algo->insert_assoc_step( poss, GeneralIterationPack::PRE_STEP, 1
			, "InitFiniteDiffReducedHessian"
			, init_red_hess_step  );
	}

	// 6/13/00: Adding steps for the computation of the exact reduced Hessian
	if( cov_.exact_reduced_hessian_ == true ) {

		// Remove the ReducedHessian steps and all of its associated pre and
		// post steps and add the exact calculation of the reduced Hessian
		Algorithm::poss_type poss;
		assert( (poss = algo->get_step_poss(ReducedHessian_name))
			!= Algorithm::DOES_NOT_EXIST );
		algo->remove_step( poss );
		algo->insert_step(
			  poss
			, ReducedHessian_name
			, new ReducedHessianExactStd_Step
		  );

		// Add the step for the exact computation of the reduced QP cross term
		assert( (poss = algo->get_step_poss(IndepDirec_name))
			!= Algorithm::DOES_NOT_EXIST );
		algo->insert_step(
			  poss
			, "ReducedQPCrossTerm"
			, new CrossTermExactStd_Step
		  );

		// Add the calculation of the dampening parameter for the cross term.
		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			DampenCrossTermStd_Step
				*zeta_step = new DampenCrossTermStd_Step;
			// ToDo: set options from stream
			algo->insert_assoc_step( poss+1, GeneralIterationPack::POST_STEP, 1
				, "DampenReducedQPCrossTerm"
				, zeta_step  );
		}

		// Change the type of the iteration quantity for rHL
		typedef GeneralIterationPack::IterQuantityAccessContinuous<
					MatrixSymPosDefCholFactor >
				iq_rHL_concrete_t;
		typedef GeneralIterationPack::IterQuantityAccessDerivedToBase< MatrixWithOp
					, MatrixSymPosDefCholFactor >
				iq_rHL_t;
		algo->state().get_iter_quant(algo->state().get_iter_quant_id(rHL_name))
			= AlgorithmState::IQ_ptr( new iq_rHL_t( new iq_rHL_concrete_t(1,rHL_name) ) );

		// That's it man!
	}

	// 4/10/01: Set the QP solver and tester objects.
	if( feasibility_step_strategy.get() ) {
		feasibility_step_strategy->set_qp_solver(qp_solver);
		feasibility_step_strategy->set_qp_tester(qp_tester);
	}

*/

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
