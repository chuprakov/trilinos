// /////////////////////////////////////////////////////////////////////////
// rSQPAlgo_ConfigMamaJama.cpp

//#define RELEASE_TRACE

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include <iostream>

#include "Misc/include/debug.h"

#include "rSQPAlgo_ConfigMamaJama.h"
#include "../../include/rSQPAlgo.h"
#include "../../include/rSQPAlgoContainer.h"
#include "../../include/rSQPStateContinuousStorage.h"
#include "../../include/rSQPStateContinuousStorageMatrixWithOpCreatorAggr.h"
#include "SparseLinAlgPack/include/COOMatrixWithPartitionedViewSubclass.h"			// Hf, Gc, Hcj
#include "ConstrainedOptimizationPack/include/IdentZeroVertConcatMatrixSubclass.h"	// Y
#include "ConstrainedOptimizationPack/include/DenseIdentVertConcatMatrixSubclass.h"	// Z
#include "ConstrainedOptimizationPack/include/ZAdjointFactMatrixSubclass.h"			// Z
#include "SparseLinAlgPack/include/COOMatrixPartitionViewSubclass.h"				// U
#include "SparseLinAlgPack/include/GenMatrixSubclass.h"								// V
#include "ConstrainedOptimizationPack/include/SymInvCholMatrixSubclass.h"			// rHL
#include "ConstrainedOptimizationPack/include/SymLBFGSMatrixSubclass.h"				// rHL
#include "ConstrainedOptimizationPack/include/SymMatrixSubclass.h"					// rHL

#include "NLPInterfacePack/include/NLPReduced.h"
#include "SparseSolverPack/include/COOBasisSystem.h"
#include "SparseSolverPack/include/MA28SparseCOOSolverCreator.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearchArmQuad_Strategy.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPL1.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPModL1.h"
#include "../../include/std/DecompositionSystemVarReductStd.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateDirect.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateAdjoint.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearchArmQuad_StrategySetOptions.h"

#include "../../include/std/ReducedQPSolverCheckOptimality.h"
#include "../../include/std/ReducedQPSolverQPOPTSOLStd.h"
#include "../../include/std/ReducedQPSolverQPSOL.h"
#include "../../include/std/ReducedQPSolverQPOPT.h"
#include "../../include/std/ReducedQPSolverQPKWIKNEWStd.h"
#include "../../include/std/ReducedQPSolverQPKWIKNEW.h"
//#include "../../include/std/QPSolverWithBoundsVE09Std.h"
//#include "../../include/std/QPSolverWithBoundsValidateInput.h"
#include "../../include/std/QPMixedFullReducedSolverCheckOptimality.h"
#include "../../include/std/QPSCPD/QPMixedFullReducedQPSCPDSolver.h"
#include "../../include/std/QPSCPD/QPMixedFullReducedQPSCPDSolverSetOptions.h"

#include "../../include/std/rSQPAlgorithmStepNames.h"

#include "../../include/std/DecompositionSystemVarReductStd.h"
#include "../../include/std/EvalNewPointStd_Step.h"
#include "../../include/std/ReducedGradientStd_Step.h"
#include "../../include/std/InitFinDiffReducedHessian_Step.h"
#include "../../include/std/InitFinDiffReducedHessian_StepSetOptions.h"
#include "../../include/std/ReducedHessianBFGSStd_Step.h"
#include "../../include/std/DepDirecStd_Step.h"
#include "../../include/std/CheckBasisFromPy_Step.h"
#include "../../include/std/IndepDirecWithoutBounds_Step.h"
#include "../../include/std/IndepDirecExact_Step.h"
#include "../../include/std/SearchDirecMixedFullReduced_Step.h"
#include "../../include/std/CalcDFromYPYZPZ_Step.h"
#include "../../include/std/CalcZpzFromDYPY_Step.h"
#include "../../include/std/LineSearchFailureNewBasisSelection_Step.h"
#include "../../include/std/NewBasisSelectionStd_Strategy.h"
#include "../../include/std/LineSearchFullStep_Step.h"
#include "../../include/std/LineSearchDirect_Step.h"
#include "../../include/std/LineSearch2ndOrderCorrect_Step.h"
#include "../../include/std/LineSearch2ndOrderCorrect_StepSetOptions.h"
#include "../../include/std/LineSearchWatchDog_Step.h"
#include "../../include/std/LineSearchFullStepAfterKIter_Step.h"
#include "../../include/std/CalcLambdaIndepStd_AddedStep.h"
#include "../../include/std/CalcReducedGradLagrangianStd_AddedStep.h"
#include "../../include/std/CheckConvergenceStd_AddedStep.h"
#include "../../include/std/CheckSkipBFGSUpdateStd_Step.h"
#include "../../include/std/CheckSkipBFGSUpdateNoPy_Step.h"
#include "../../include/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.h"
#include "../../include/std/MeritFunc_PenaltyParamUpdateWithMult_AddedStep.h"
#include "../../include/std/MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.h"
#include "../../include/std/MeritFunc_ModifiedL1LargerSteps_AddedStep.h"
#include "../../include/std/MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.h"
#include "../../include/std/ActSetStats_AddedStep.h"
#include "../../include/std/NumFixedDepIndep_AddedStep.h"
#include "../../include/std/PrintReducedQPError_AddedStep.h"

#include "../../include/std/act_set_stats.h"
#include "../../include/std/qp_solver_stats.h"
#include "../../include/std/quasi_newton_stats.h"

#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateDirect.h"
#include "SparseLinAlgPack/include/sparse_bounds.h"

namespace ReducedSpaceSQPPack {

rSQPAlgo_ConfigMamaJama::rSQPAlgo_ConfigMamaJama(
		  ReferenceCountingPack::ref_count_ptr<BasisSystem> basis_sys_ptr
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			Gc_iq_creator_ptr
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			U_iq_creator_ptr	)
	: algo_(0)
		, algo_cntr_(0)
		, in_destructor_(false)
		, basis_sys_ptr_(basis_sys_ptr)
		, Gc_iq_creator_ptr_(Gc_iq_creator_ptr)
		, U_iq_creator_ptr_(U_iq_creator_ptr)
		, qp_solver_type_(QPOPT)
		, factorization_type_(AUTO_FACT)
		, bigM_(-1.0)
		, max_basis_cond_change_frac_(-1.0)
		, warm_start_frac_(-1.0)
		, merit_function_type_(MERIT_FUNC_MOD_L1_INCR)
		, line_search_method_(LINE_SEARCH_DIRECT)
		, use_line_search_correct_kkt_tol_(-1.0)
		, full_steps_after_k_(-1)
		, quasi_newton_(QN_AUTO)
		, hessian_initialization_(INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS)
		, quasi_newton_dampening_(false)
		, num_lbfgs_updates_stored_(7)
		, lbfgs_auto_scaling_(true)
		, max_dof_quasi_newton_dense_(1000)
		, check_results_(false)
{}

rSQPAlgo_ConfigMamaJama::~rSQPAlgo_ConfigMamaJama() {

	in_destructor_ = true;

	if(algo_) {
		// If the config object is deleted before the algo container (algo_cntr)
		// before algo_cntr has been deconfigured then
		// either #this# was allocatted on the stack after algo_cntr in which case
		// we don't want algo_cntr call or delete #this# in its destructor, or
		// #this# was dynamically allocated and is now being deleted by someone
		// else by mistake.

		// The expression (algo_ && algo_cntr_) can not be true
		assert( algo_cntr_ );
		// If algo_cntr_ owns #this# then this->decomfig_algo_cntr()
		// should have been called first so you should not be here.		
		assert( !algo_cntr_->get_config().has_ownership() );

		// Must defig the algo and algo_cntr because I will not be around to call
		// after this.
		algo_cntr_->set_algo(0);
		algo_cntr_->set_config(0);
		delete algo_;
	}
}

void rSQPAlgo_ConfigMamaJama::warm_start_frac(value_type warm_start_frac) {
	if( warm_start_frac <= 0.0 || warm_start_frac > 1.0 )
		throw std::logic_error( "rSQPAlgo_ConfigMamaJama::warm_start_frac(frac) : "
			"Error, frac must be in the range (0,1]"								);
	warm_start_frac_ = warm_start_frac;
}

void rSQPAlgo_ConfigMamaJama::config_algo_cntr(rSQPAlgoContainer& algo_cntr) {

//	TRACE0( "\n*** rSQPAlgo_ConfigMamaJama::config_algo_cntr(algo_cntr) called\n" );
//	TRACE0( "\n*** Before algo is configured\n" );
//	print_state();

	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	
	if(algo_cntr_) {
		if(algo_cntr_ != &algo_cntr)
			throw std::logic_error("rSQPAlgo_ConfigMamaJama::config_algo_cntr(): rSQPAlgo_ConfigMamaJama can only configure one rSQPAlgoContainer and I have already configured mine");
		else
			return; // same one so you are done
	}

	// ////////////////////////////////////////////////////////////
	// A. Remember the algorithm container

	// remember who the algo container is
	algo_cntr_ = &algo_cntr;

	// /////////////////////////////////////////////////////////////////////////
	// B. Create an algo object, give to algo_cntr, then give algo_cntr to algo

#ifdef RELEASE_TRACE
	std::cout << "\n*** Create the algo object ...\n";
#endif

	assert(algo_ = new rSQPAlgo);

#ifdef RELEASE_TRACE
	std::cout << "\n*** algo_cntr_->set_algo(algo_) ...\n";
#endif

	algo_cntr_->set_algo(algo_);

#ifdef RELEASE_TRACE
	std::cout << "\n*** algo_->set_algo_cntr(algo_cntr_) ...\n";
#endif

	algo_->set_algo_cntr(algo_cntr_);

	// /////////////////////////////////////////////
	// C. Configure algo

	// /////////////////////////////////////////////////////
	// C.0 Set the nlp and track objects

#ifdef RELEASE_TRACE
	std::cout << "\n*** Set the NLP and track objects ...\n";
#endif

	algo_->set_nlp( algo_cntr_->get_nlp().get() );
	algo_->set_track( rcp::rcp_implicit_cast<AlgorithmTrack>(algo_cntr_->get_track()) );

	// Determine whether to use direct or adjoint factorization
	switch(factorization_type_) {
		case DIRECT_FACT:
			fact_type_ = DIRECT_FACT;
			break;
		case ADJOINT_FACT:
			fact_type_ = ADJOINT_FACT;
			break;
		case AUTO_FACT:
		{
			if( ! algo_->nlp().has_bounds()  || qp_solver_type_ == VE09
					|| qp_solver_type_ == QPSCPD )
			{
				fact_type_ = ADJOINT_FACT;
			}
			else {
				// If the total number of bounds is less than the number
				// of degrees of freedom then the adjoint factorization
				// will be faster.
				NLPReduced &nlp = algo_->nlp();
				if(!nlp.is_initialized()) nlp.initialize();
				using SparseLinAlgPack::num_bounds;
				if( num_bounds( nlp.xl(), nlp.xu() ) < nlp.n() - nlp.r() )
					fact_type_ = ADJOINT_FACT;
				else
					fact_type_ = DIRECT_FACT;
			}
			break;
		}
	}

	// /////////////////////////////////////////////////////
	// C.1. Create and set the state object

#ifdef RELEASE_TRACE
	std::cout << "\n*** Create and set the state object ...\n";
#endif

	{
		// set the matrix creator (factory)
		typedef rSQPStateContinuousStorageMatrixWithOpCreator
			matrix_creator_t;
		typedef rSQPStateContinuousStorageMatrixWithOpCreatorAggr
			matrix_creator_aggr_t;
		typedef rSQPStateContinuousStorageMatrixWithOpCreatorAggr::matrix_iqa_creator_ptr_t
			matrix_iqa_creator_ptr_t;

		matrix_iqa_creator_ptr_t matrix_iqa_creators[matrix_creator_t::num_MatrixWithOp] =
		{
			new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>(),						// Q_Hf
			Gc_iq_creator_ptr_.get()
				? Gc_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<COOMatrixWithPartitionedViewSubclass>(),		// Q_Gc
			new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>(),						// Q_Hcj
			new IterQuantMatrixWithOpCreatorContinuous<IdentZeroVertConcatMatrixSubclass>(),		// Q_Y
			0,																						// Q_Z
			U_iq_creator_ptr_.get()
				? U_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<COOMatrixPartitionViewSubclass>(),			// Q_U
			new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>(),						// Q_V
			0																						// Q_rHL
		};

		// Q_Z
		switch(fact_type_) {
			case DIRECT_FACT:
				matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_Z]
					= new IterQuantMatrixWithOpCreatorContinuous<DenseIdentVertConcatMatrixSubclass>();
				break;
			case ADJOINT_FACT:
				matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_Z]
					= new IterQuantMatrixWithOpCreatorContinuous<ZAdjointFactMatrixSubclass>();
				break;
			default:
				assert(0);
		}
		
		// Q_rHL

		// Decide if we could use limited memory BFGS or not?
		bool use_limited_memory = false;
		const size_type
			n = algo_->nlp().n(),
			m = algo_->nlp().m();
		switch( quasi_newton_ ) {
			case QN_AUTO: {
				if( n - m > max_dof_quasi_newton_dense_ )
					use_limited_memory = true;
				break;
			}
			case QN_BFGS:
				use_limited_memory = false;
				break;
			case QN_LBFGS:
				use_limited_memory = true;
				break;
			default:
				assert(0);	// only local programming error
		}

		if( ! algo_->nlp().has_bounds() ) {
			// Here we just need to solve for linear systems with rHL
			// and we can use a limited memory method or update the dense
			// factors directly.
			matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
				= use_limited_memory
					? static_cast<IterQuantMatrixWithOpCreator*>(
						new IterQuantMatrixWithOpCreatorContinuous<
							SymLBFGSMatrixSubclass>() )
					: static_cast<IterQuantMatrixWithOpCreator*>( 
						new IterQuantMatrixWithOpCreatorContinuous<
							SymInvCholMatrixSubclass>() );
		}
		else {
			switch(qp_solver_type_) {
				case QPSOL:
				case QPOPT:
				case VE09:
					matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
						= use_limited_memory
							? static_cast<IterQuantMatrixWithOpCreator*>(
								new IterQuantMatrixWithOpCreatorContinuous<
									SymLBFGSMatrixSubclass>() )
							: static_cast<IterQuantMatrixWithOpCreator*>( 
								new IterQuantMatrixWithOpCreatorContinuous<
									SymMatrixSubclass>() );
					break;
				case QPKWIK:
					matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
						= new IterQuantMatrixWithOpCreatorContinuous<
									SymInvCholMatrixSubclass>();
					use_limited_memory = false;
					break;
				case QPSCPD:
					matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
						= use_limited_memory
							? static_cast<IterQuantMatrixWithOpCreator*>(
								new IterQuantMatrixWithOpCreatorContinuous<
									SymLBFGSMatrixSubclass>() )
							: static_cast<IterQuantMatrixWithOpCreator*>( 
								new IterQuantMatrixWithOpCreatorContinuous<
									SymInvCholMatrixSubclass>() );
					break;
				default:
					assert(0);
			}
		}
				
		// set the number of storage locations
		typedef rSQPStateContinuousStorage state_t;
		size_type storage_num[state_t::num_quantities];
		state_t::set_default_storage_num(storage_num);

		storage_num[state_t::Q_x]				= 2;
		storage_num[state_t::Q_f]				= 2;
		storage_num[state_t::Q_c]				= 2;
		storage_num[state_t::Q_norm_inf_c]		= 2;

		storage_num[state_t::Q_rGL]				= 2;
		storage_num[state_t::Q_norm_2_rGL]		= 2;
		storage_num[state_t::Q_norm_inf_rGL]	= 2;
		storage_num[state_t::Q_lambda]			= 2;

		storage_num[state_t::Q_d]				= 2;
		storage_num[state_t::Q_norm_2_d]		= 2;
		storage_num[state_t::Q_norm_inf_d]		= 2;

		storage_num[state_t::Q_rGf]				= 2;

		storage_num[state_t::Q_alpha]			= 2;
		storage_num[state_t::Q_mu]				= 2;
		storage_num[state_t::Q_phi]				= 2;
		storage_num[state_t::Q_nu]				= 2;

		// create a new state object
		algo_->set_state( new state_t( new matrix_creator_aggr_t(matrix_iqa_creators), storage_num ) );

		// Now set the number of LBFGS update vectors
		if(use_limited_memory) {
			// Here we will set the number of update vectors to store by setting
			// it for the k-1 iteration.  That way the k iteration will still not
			// be updated and therefore the regular initializations will still
			// be performed.  This is a little bit of a hack but it should work.
			SymLBFGSMatrixSubclass *_rHL
				= dynamic_cast<SymLBFGSMatrixSubclass*>(&algo_->rsqp_state().rHL().set_k(-1));
			if(_rHL) {
				_rHL->set_num_updates_stored( num_lbfgs_updates_stored_ );
				_rHL->auto_rescaling( lbfgs_auto_scaling_ );
				algo_->rsqp_state().rHL().set_not_updated(-1);
			}
		}
	}

	// /////////////////////////////////////////////////////
	// C.2. Create and set the decomposition system object

#ifdef RELEASE_TRACE
	std::cout << "\n*** Create and set the decomposition system object ...\n";
#endif

	// ToDo: guard from an exception being thrown from these
	// constructors.

	if( !basis_sys_ptr_.get() ) {

		MA28SparseCOOSolverCreator*
			solver_creator = new MA28SparseCOOSolverCreator(
					new SparseSolverPack::MA28SparseCOOSolverSetOptions
					, const_cast<OptionsFromStreamPack::OptionsFromStream*>(options_)
				);
		
		basis_sys_ptr_ = new COOBasisSystem(solver_creator, true);

	}

	DecompositionSystemVarReductImpNode* decomp_sys_aggr = 0;
	
	switch(fact_type_) {
		case DIRECT_FACT:
			decomp_sys_aggr = new DecompositionSystemCoordinateDirect(
									basis_sys_ptr_.release(), true);
			break;
		case ADJOINT_FACT:
			decomp_sys_aggr = new DecompositionSystemCoordinateAdjoint(
									basis_sys_ptr_.release(), true);
			break;
	}
	
	DecompositionSystemVarReductStd*
		decomp_sys = new DecompositionSystemVarReductStd( false, algo_, decomp_sys_aggr );

	algo_->rsqp_state().set_decomp_sys(decomp_sys);

	// /////////////////////////////////////////////////////
	// C.3  Create and set the step objects

#ifdef RELEASE_TRACE
	std::cout << "\n*** Create and set the step objects ...\n";
#endif

	{
		typedef ref_count_ptr<MeritFuncNLP> merit_func_ptr_t;
		merit_func_ptr_t merit_func;

		int step_num = 0;

		// Set the standard set of step objects.
		// This default algorithm is for NLPs with no variable bounds.

		// (1) EvalNewPoint
		algo_->insert_step( ++step_num, EvalNewPoint_name, new EvalNewPointStd_Step );

		// (2) ReducedGradient
		algo_->insert_step( ++step_num, ReducedGradient_name, new ReducedGradientStd_Step );

		// (3) Calculate Reduced Gradient of the Lagrangian
		algo_->insert_step( ++step_num, CalcReducedGradLagrangian_name, new CalcReducedGradLagrangianStd_AddedStep );

		// (4)	Calculate the Lagrange multipliers for the independent constraints.
		// 		These are computed here just in case the algorithm converges and we need to
		// 		report these multipliers to the NLP.
		algo_->insert_step( ++step_num, CalcLambdaIndep_name, new  CalcLambdaIndepStd_AddedStep		);

		// (5) Check for convergence
		algo_->insert_step( ++step_num, CheckConvergence_name, new CheckConvergenceStd_AddedStep );

		// (6) ReducedHessian
		ref_count_ptr<ReducedHessianBFGSStd_Step>
			bfgs_updater = new ReducedHessianBFGSStd_Step(quasi_newton_dampening_);
		algo_->insert_step( ++step_num, ReducedHessian_name
			, rcp::rcp_implicit_cast<AlgorithmStep>(bfgs_updater) );

		// (6.-1) CheckSkipBFGSUpdate
		algo_->insert_assoc_step(
			  step_num
			, GeneralIterationPack::PRE_STEP
			, 1
			, CheckSkipBFGSUpdate_name
			, new CheckSkipBFGSUpdateStd_Step
		  );

		// (7) DepDirec
		algo_->insert_step( ++step_num, DepDirec_name, new  DepDirecStd_Step );

		// (8) IndepDirec
		algo_->insert_step( ++step_num, IndepDirec_name, new  IndepDirecWithoutBounds_Step );

		// (9) CalcDFromYPYZPZ

		algo_->insert_step( ++step_num, "CalcDFromYPYZPZ", new CalcDFromYPYZPZ_Step );

		// (10) LineSearch

		if( line_search_method_ == LINE_SEARCH_NONE ) {

			algo_->insert_step(
				  ++step_num
				, LineSearch_name
				, new LineSearchFullStep_Step
			  );
			
		}
		else {

			using ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy;
			using ConstrainedOptimizationPack::MeritFuncNLPL1;
			using ConstrainedOptimizationPack::MeritFuncNLPModL1;

			switch( merit_function_type_ ) {
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
			ls_options_setter.set_options( *options_ );
			LineSearch_Step
				*line_search_step = 0;

			switch( line_search_method_ ) {
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
					switch( algo_cntr_->iteration_info_output() ) {
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
					// Create the line search object for the newton iterations
					// and set its options from the option from stream object.
					DirectLineSearchArmQuad_Strategy
						*direct_ls_newton = new  DirectLineSearchArmQuad_Strategy(
								5	// max_iter
							);
					ConstrainedOptimizationPack::DirectLineSearchArmQuad_StrategySetOptions
						direct_ls_opt_setter( direct_ls_newton
							, "DirectLineSearchArmQuad2ndOrderCorrectNewton" );
					direct_ls_opt_setter.set_options( *options_ );
					// Create the step object and set its options from the options object.
					LineSearch2ndOrderCorrect_Step
						*_line_search_step = new LineSearch2ndOrderCorrect_Step(
							  direct_line_search
							, rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
							, direct_ls_newton
							, direct_line_search->eta()
							, newton_olevel
						);
					LineSearch2ndOrderCorrect_StepSetOptions
						ls_opt_setter( _line_search_step );
					ls_opt_setter.set_options( *options_ );
					//
					line_search_step = _line_search_step;
					break;
				}
				case LINE_SEARCH_WATCHDOG: {
					LineSearchWatchDog_Step
						*_line_search_step = new LineSearchWatchDog_Step(
							  direct_line_search
							, rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
							, use_line_search_correct_kkt_tol_ >= 0.0
								? use_line_search_correct_kkt_tol_ : 0.0
							, direct_line_search->eta()
						);
					line_search_step = _line_search_step;
					break;
				}
			}

			algo_->insert_step(
				  ++step_num
				, LineSearch_name
				, new LineSearchFailureNewBasisSelection_Step( 
					  line_search_step
					, new NewBasisSelectionStd_Strategy(decomp_sys)
				  )
			  );

			Algorithm::poss_type
				pre_step_i = 0;

			// (10.-?) Update the penalty parameter for the Merit function.
			MeritFunc_PenaltyParamUpdate_AddedStep
				*param_update_step = 0;

			switch( merit_function_type_ ) {
				case MERIT_FUNC_L1:
					param_update_step = new  MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
											rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
					break;
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
			ppu_options_setter.set_options( *options_ );

			algo_->insert_assoc_step(	  step_num
										, GeneralIterationPack::PRE_STEP
										, ++pre_step_i
										, "MeritFunc_PenaltyParamUpdate"
										, param_update_step // give control over memory!
									);
			
			// (10.-?)	Compute the full step before the linesearch
			algo_->insert_assoc_step(	  step_num
										, GeneralIterationPack::PRE_STEP
										, ++pre_step_i
										, "LineSearchFullStep"
										, new  LineSearchFullStep_Step		);

			// (10.-?) Increase the penalty parameters to get a larger step.
			if( merit_function_type_ == MERIT_FUNC_MOD_L1_INCR ) {

				MeritFunc_ModifiedL1LargerSteps_AddedStep
					*_added_step = new MeritFunc_ModifiedL1LargerSteps_AddedStep(
						   rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func)
						 , direct_line_search->eta() );

				MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions
					options_setter( _added_step );
				options_setter.set_options( *options_ );

				algo_->insert_assoc_step(	  step_num
											, GeneralIterationPack::PRE_STEP
											, ++pre_step_i
											, "MeritFunc_ModifiedL1LargerSteps"
											, _added_step // give control over memory!
										);
			}
		}
	
#ifdef RELEASE_TRACE
		std::cout << "\n*** Reconfigure steps for NLP with bounds ...\n";
#endif

		if( algo_cntr_->nlp().has_bounds() ) {
			
			// If the NLP has bounds on the variables then we must replace
			// IndepDirec and move the convergence check to the end (along
			// with the calculation of the reduced gradient of the lagrangian).
			
			typedef ref_count_ptr<ReducedQPSolver> qp_solver_ptr_t;
			qp_solver_ptr_t qp_solver;

			switch(qp_solver_type_) {
				case QPOPT:
				case QPSOL:
				{
					ReducedQPSolverQPOPTSOL*  _qp_solver = 0;
					if(qp_solver_type_ == QPOPT)
						_qp_solver = new ReducedQPSolverQPOPT;
					else
						_qp_solver = new ReducedQPSolverQPSOL;
					
					if( bigM_ > 0.0 )
						_qp_solver->bigM(bigM_);	

					// If we are going to dump all of the iteration quantites we might as
					// well just print out the mappings to QPOPT and QPSOL.
					if(algo_cntr_->iteration_info_output() == PRINT_ITERATION_QUANTITIES)
						mapped_qp_file_ = mapped_qp_file_ptr_t( new std::ofstream("mapped_qp.txt") );
					_qp_solver->qp_mapping_output( mapped_qp_file_.get() );
					qp_solver = qp_solver_ptr_t( new ReducedQPSolverQPOPTSOLStd(_qp_solver,algo_) );
					break;
				}
				case QPKWIK:
				{
					ReducedQPSolverQPKWIKNEW*
						_qp_solver = new ReducedQPSolverQPKWIKNEW;
					if(algo_cntr_->iteration_info_output() == PRINT_ITERATION_QUANTITIES)
						_qp_solver->create_qpkwiknew_file(true);
					ReducedQPSolverQPKWIKNEWStd*
						__qp_solver = new ReducedQPSolverQPKWIKNEWStd(_qp_solver,algo_);
					if( warm_start_frac_ > 0.0 )
						__qp_solver->warm_start_frac(warm_start_frac_);
					qp_solver = qp_solver_ptr_t( __qp_solver );
					break;
				}
			}

			Algorithm::poss_type poss;

			poss = algo_->get_step_poss(IndepDirec_name);
			algo_->replace_step(
						poss
						,new  IndepDirecExact_Step( 
							new ReducedQPSolverCheckOptimality(
								qp_solver, algo_, check_results_ )
						  )
				    );

			Algorithm::step_ptr_t calc_rgrad_lagr, calc_lambda, check_conv;

			// Remove and save CalcReducedGradLagr..., CalcLambdaIndep... and CheckConv...
			



			poss			= algo_->get_step_poss(CalcReducedGradLagrangian_name);
			calc_rgrad_lagr	= algo_->get_step(poss);
			algo_->remove_step(poss);

			poss			= algo_->get_step_poss(CalcLambdaIndep_name);
			calc_lambda		= algo_->get_step(poss);
			algo_->remove_step(poss);

			poss			= algo_->get_step_poss(CheckConvergence_name);
			check_conv		= algo_->get_step(poss);
			algo_->remove_step(poss);

			// Add them before LineSearch
	
			poss		= algo_->get_step_poss(LineSearch_name);
			algo_->insert_step( poss++, CalcReducedGradLagrangian_name, calc_rgrad_lagr );
			algo_->insert_step( poss++, CalcLambdaIndep_name, calc_lambda );
			algo_->insert_step( poss++, CheckConvergence_name, check_conv );

		}

		// Added to incorrperate VE09.
		// Latter modified to also incorperate QPSCPD
		if( qp_solver_type_ == VE09 || qp_solver_type_ == QPSCPD ) {

			// Here I will replace the IndepDirec step object and insert
			// the step object that will compute the full d and lambda_indep.
			// Then I will insert a step to compute Zpz = d - Ypy so that
			// I can continue to use the BFGS updating checker.
			// I can also remove the step that computes lambda_indep
			// since it is already computed here.

			// note: ref_count_ptr<> manages memory here.
			using QPSCPDPack::QPMixedFullReducedQPSCPDSolver;
			QPMixedFullReducedQPSCPDSolver
				*qp_solver = new QPMixedFullReducedQPSCPDSolver(algo_);
			qp_solver->set_pick_violated_policy(
					ReducedSpaceSQPPack::QPSCPDPack::ConstraintsVarBoundsRelaxed::MOST_VIOLATED );
			if( bigM_ > 0.0 )
				qp_solver->set_M(bigM_);
			if( warm_start_frac_ > 0.0 )
				qp_solver->set_warm_start_frac(warm_start_frac_);
			using ConstrainedOptimizationPack::QPSCPD;	
			QPSCPD::EOutputLevel	qp_olevel;
			switch( algo_cntr_->iteration_info_output() ) {
				case PRINT_NOTHING:
					qp_olevel = QPSCPD::NO_OUTPUT;
					break;
				case PRINT_BASIC_ALGORITHM_INFO:
					qp_olevel = QPSCPD::OUTPUT_BASIC_INFO;
					break;
				case PRINT_ALGORITHM_STEPS:
					qp_olevel = QPSCPD::OUTPUT_BASIC_INFO;
					break;
				case PRINT_ACTIVE_SET:
					qp_olevel = QPSCPD::OUTPUT_ITER_SUMMARY;
					break;
				case PRINT_VECTORS:
					qp_olevel = QPSCPD::OUTPUT_ITER_SUMMARY;
					break;
				case PRINT_ITERATION_QUANTITIES:
					qp_olevel = QPSCPD::OUTPUT_ITER_QUANTITIES;
					break;
			}
			qp_solver->set_print_level( qp_olevel );
			if( options_ ) {
				ReducedSpaceSQPPack::QPSCPDPack::QPMixedFullReducedQPSCPDSolverSetOptions
					option_setter( qp_solver );
				option_setter.set_options( *options_ );
			}

			QPMixedFullReducedSolverCheckOptimality
				*qp_opt_checker = new QPMixedFullReducedSolverCheckOptimality(
					 qp_solver, algo_, check_results_ );

			Algorithm::poss_type poss;
			poss = algo_->get_step_poss(IndepDirec_name);

			algo_->remove_step(poss);

			algo_->insert_step(
				poss
				, SearchDirec_name
				, new SearchDirecMixedFullReduced_Step(
					qp_opt_checker
				  )
			  );

			algo_->insert_step(
				poss + 1
				, "CalcZpzFromDYPY"
				, new CalcZpzFromDYPY_Step
			  );

			poss = algo_->get_step_poss(CalcLambdaIndep_name);
			algo_->remove_step(poss);

/* VE09 is disabled for now.  Latter I will reintegrate it with what is above.

			// Set the options for VE09 QP solver
			QPSolverWithBoundsVE09 *qp_solver = new QPSolverWithBoundsVE09;
			qp_solver->create_ve09_err_mon_files(
				algo_cntr_->iteration_info_output() >= PRINT_ALGORITHM_STEPS );

			Algorithm::poss_type poss;
			poss = algo_->get_step_poss( IndepDirec_name );
			algo_->replace_step(
				poss
				, new SearchDirecMixedFullReduced_Step(
					new QPSolverWithBoundsValidateInput(
						new QPSolverWithBoundsVE09Std(
							qp_solver
							, algo_
						)
					)
				  )
			  );

*/

/* This algorithm has been disabled because it did not work (I think).
			// Here we must remove DepDirec and IndepDirec and replace them with
			// SearchDirec.  Since the multipliers for the constraints will
			// also be calculated in the calculation of the full
			// space search directioin we can remove CalcLambdaIndep.
			// We also need to replace how the BFGS update is checked since
			// py is no longer calculated.

			// Replace calculation of search direction

			algo_->remove_step( algo_->get_step_poss( DepDirec_name ) );
			algo_->remove_step( algo_->get_step_poss( IndepDirec_name ) );
			algo_->insert_step(
				algo_->get_step_poss( ReducedHessian_name ) + 1	// after reduced hessian calc.
				, SearchDirec_name
				, new SearchDirecMixedFullReduced_Step(
					new QPSolverWithBoundsValidateInput(
						new QPSolverWithBoundsVE09Std(
							new QPSolverWithBoundsVE09
							, algo_
						)
					)
				  )
			  );

			// Remove the CalcLambdaIndep step

			Algorithm::poss_type poss;
			poss = algo_->get_step_poss( LineSearch_name );
			algo_->remove_assoc_step(
				poss
				, GeneralIterationPack::PRE_STEP
				, algo_->get_assoc_step_poss(
					poss
					, GeneralIterationPack::PRE_STEP
					, CalcLambdaIndep_name
				  )
			  );

			// Replace Check for skipping BFGS update

			poss = algo_->get_step_poss( ReducedHessian_name );

			Algorithm::poss_type assoc_poss;
			
			algo_->remove_assoc_step(
				poss
				, GeneralIterationPack::PRE_STEP
				, assoc_poss = algo_->get_assoc_step_poss(
					poss
					, GeneralIterationPack::PRE_STEP
					, CheckSkipBFGSUpdate_name			  )
			  );

			algo_->insert_assoc_step(
				poss
				, GeneralIterationPack::PRE_STEP
				, assoc_poss
				, CheckSkipBFGSUpdate_name
				, new CheckSkipBFGSUpdateNoPy_Step
			  );
*/

		}
	}

	// 8/30/99: Add the iteration quantities for the QPSolverStats and 
	// ActSetStats and the added step that will calculate the
	// active set updates.

	{
		// Add active set and QP statistics to state object
		algo_->state().set_iter_quant( act_set_stats_name
			, new IterQuantityAccessContinuous<ActSetStats>( 1, act_set_stats_name ) );
		algo_->state().set_iter_quant( qp_solver_stats_name
			, new IterQuantityAccessContinuous<QPSolverStats>( 1, qp_solver_stats_name ) );

		// Insert active set computational step into algorithm
		Algorithm::poss_type
			poss = algo_->get_step_poss( IndepDirec_name );
		if( poss == Algorithm::DOES_NOT_EXIST ) {
			// If we are not using IndepDirec to compute the active set then
			// we must be using SearchDirec
			assert( poss = algo_->get_step_poss( SearchDirec_name ) );		
		}
		Algorithm::poss_type
			nas = algo_->num_assoc_steps( poss, GeneralIterationPack::POST_STEP );
		algo_->insert_assoc_step( poss, GeneralIterationPack::POST_STEP, nas+1
			, "ActiveSetStatistics", new ActSetStats_AddedStep );
		
		// Output the number of fixed depenent and indepent variables.
		algo_->insert_assoc_step( poss, GeneralIterationPack::POST_STEP, nas+2
			, "CountNumFixedDepIndep", new NumFixedDepIndep_AddedStep );

		// Output of KKT conditions of reduced QP subproblem
		algo_->insert_assoc_step( poss, GeneralIterationPack::POST_STEP, nas+3
			, "PrintReducedQPError", new PrintReducedQPError_AddedStep( print_qp_error_ ) );

	}

	// 10/15/99: Add basis checking step
	{
		CheckBasisFromPy_Step
			*basis_check_step = new CheckBasisFromPy_Step(
				new NewBasisSelectionStd_Strategy(decomp_sys) );
		if( max_basis_cond_change_frac_ > 0.0 )
			basis_check_step->max_basis_cond_change_frac( max_basis_cond_change_frac_ );
		Algorithm::poss_type poss;
		assert(poss = algo_->get_step_poss( DepDirec_name ) );
		algo_->insert_step(
			  poss+1
			, "CheckBasis"
			, basis_check_step
		  );
	}

	// 10/19/99: Add quasi-newton stats
	algo_->state().set_iter_quant( quasi_newton_stats_name
		, new IterQuantityAccessContinuous<QuasiNewtonStats>( 1, quasi_newton_stats_name ) );

	// 11/12/99: Add the option of using full steps after a specified number of iterations
	// if this option has been set and if we are using a linesearch method.
	if( full_steps_after_k_ > 0 && line_search_method_ != LINE_SEARCH_NONE ) {
		Algorithm::poss_type poss;
		assert(poss = algo_->get_step_poss( LineSearch_name ) );
		Algorithm::step_ptr_t
			existing_step_ptr = algo_->get_step( poss );
		// Check that the existing step is indeed a linesearch step and that
		// we can safely use static cast cast up in types.  If multiple
		// inheritance was used then this test will also fail.
		assert( dynamic_cast<LineSearch_Step*>(existing_step_ptr.get())
				==  existing_step_ptr.get() );
		rcp::ref_count_ptr<LineSearch_Step>
			existing_ls_step_ptr = rcp::rcp_static_cast<LineSearch_Step>(existing_step_ptr);
		algo_->replace_step(
			  poss
			, new LineSearchFullStepAfterKIter_Step(
				  existing_ls_step_ptr
				, full_steps_after_k_
			  )
		  );

	}

	// 12/3/99: Adding finite difference initializaiton of the reduced hessian
	if( hessian_initialization_ != INIT_HESS_IDENTITY ) {
		Algorithm::poss_type poss;
		assert(poss = algo_->get_step_poss( ReducedHessian_name ) );
		InitFinDiffReducedHessian_Step::EInitializationMethod
			init_hess;
		switch( hessian_initialization_ ) {
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
		opt_setter.set_options( *options_ );
		algo_->insert_assoc_step( poss, GeneralIterationPack::PRE_STEP, 1
			, "InitFiniteDiffReducedHessian"
			, init_red_hess_step  );
	}

	// Desprite debugging stuff
//	assert(algo_);
//	TRACE0( "\n*** After algo is configured\n" );
//	print_state();
}

void rSQPAlgo_ConfigMamaJama::deconfig_algo_cntr(rSQPAlgoContainer& algo_cntr) {

//	TRACE0( "\n*** rSQPAlgo_ConfigMamaJama::deconfig_algo_cntr(algo_cntr) called\n" );
//	TRACE0( "\n*** Before algo is deconfigured\n" );
//	print_state();

	if(in_destructor_) return;	// This funciton was called from algo_cntr.set_config()
								// which was called from ~rSQPAlgo_ConfigMamaJama()
	if( algo_cntr.get_algo() != algo_ )
		throw std::logic_error("rSQPAlgo_ConfigMamaJama::deconfig_algo_cntr(...): This is not my algo object");
	delete algo_;
	algo_cntr.set_algo(algo_ = 0);

//	TRACE0( "\n*** After algo is deconfigured\n" );
//	print_state();

}

void rSQPAlgo_ConfigMamaJama::init_algo(rSQPAlgoInterface& _algo)
{

//	TRACE0( "\n*** rSQPAlgo_ConfigMamaJama::init_algo(_algo) called\n" );
//	TRACE0( "\n*** Before algo is initialized\n" );
//	print_state();
		
	namespace rcp = ReferenceCountingPack;

	assert( algo_ );			// The algo object should have been set

	assert( &_algo == algo_  );	// They should be the same

	rSQPAlgo	&algo	= *algo_;
	rSQPState	&state	= algo.rsqp_state();

	algo_->max_iter( algo_cntr_->max_iter() );
	algo_->rsqp_state().iteration_info_output( algo_cntr_->iteration_info_output() );
	algo_->rsqp_state().check_results( check_results_ );

	NLPReduced	&nlp	= algo.nlp();

	// (1) Initialize the nlp
	nlp.initialize();

	// (2) Set the initial x to the initial guess.
	state.x().set_k(0) = nlp.xinit();

	// (3) Get the initial basis from the nlp.  The nlp will have a basis
	// set whether or not nlp.nlp_selects_basis() returns true (See NLPReduced).
	size_type rank;
	nlp.get_basis( &state.var_perm_new(), &state.con_perm_new(), &rank );
	size_type	n = nlp.n(),
				m = nlp.m(),
				r = nlp.r();
	state.var_dep()		= Range1D(1, r);
	state.var_indep()	= Range1D(r + 1, n);
	state.con_indep()	= Range1D(1, r);
	state.con_dep()		= (r < m) ? Range1D(r + 1, m) : Range1D(m + 1, m + 1);

	state.initialize_fast_access();

	// set the first step
	algo.do_step_first(1);

	// The rest of the algorithm will initialize itself
}

void rSQPAlgo_ConfigMamaJama::print_algorithm(std::ostream& out) const {
	out
		<< "\n"
		<< "*************************************\n"
		<< "*** rSQPAlgo_ConfigMamaJama setup ***\n"
		<< "*************************************\n";
	
	assert(algo_);

	algo_->print_algorithm(out);
}

void rSQPAlgo_ConfigMamaJama::print_state() const {
//	TRACE0("\n*** rSQPAlgo_ConfigMamaJama State ***\n" );
//	TRACE1("algo_ = 0x%x\n", algo_ );
//	TRACE1("algo_cntr_ = 0x%x\n", algo_cntr_ );
}

}	// end namespace ReducedSpaceSQPPack 