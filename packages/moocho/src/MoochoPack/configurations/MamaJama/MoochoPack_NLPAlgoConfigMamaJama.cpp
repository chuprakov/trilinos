// /////////////////////////////////////////////////////////////////////////
// rSQPAlgo_ConfigMamaJama.cpp

//#define RELEASE_TRACE

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include <sstream>
#include <typeinfo>

#include <iostream>

#include "Misc/include/debug.h"

#include "rSQPAlgo_ConfigMamaJama.h"
#include "../../include/rSQPAlgo.h"
#include "../../include/rSQPAlgoContainer.h"
#include "../../include/rSQPStateContinuousStorage.h"
#include "../../include/rSQPStateContinuousStorageMatrixWithOpCreatorAggr.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/COOMatrixWithPartitionedViewSubclass.h"			// HL, Gc, Hcj
#include "ConstrainedOptimizationPack/include/IdentZeroVertConcatMatrixSubclass.h"	// Y
#include "ConstrainedOptimizationPack/include/DenseIdentVertConcatMatrixSubclass.h"	// Z
#include "ConstrainedOptimizationPack/include/ZAdjointFactMatrixSubclass.h"			// Z
#include "SparseLinAlgPack/include/COOMatrixPartitionViewSubclass.h"				// U
#include "SparseLinAlgPack/include/GenMatrixSubclass.h"								// V
#include "ConstrainedOptimizationPack/include/SymInvCholMatrixSubclass.h"			// rHL
#include "ConstrainedOptimizationPack/include/SymLBFGSMatrixSubclass.h"				// rHL
#include "ConstrainedOptimizationPack/include/SymMatrixSubclass.h"					// rHL

#include "ReducedSpaceSQPPack/include/NLPrSQPTailoredApproach.h"
#include "ReducedSpaceSQPPack/include/NLPrSQPTailoredApproachTester.h"
#include "ReducedSpaceSQPPack/include/NLPrSQPTailoredApproachTesterSetOptions.h"

#include "NLPInterfacePack/test/NLPFirstDerivativesTester.h"
#include "NLPInterfacePack/test/NLPFirstDerivativesTesterSetOptions.h"

#include "NLPInterfacePack/include/NLPReduced.h"
#include "SparseSolverPack/include/COOBasisSystem.h"
#include "SparseSolverPack/include/MA28SparseCOOSolverCreator.h"
#include "SparseSolverPack/include/MA48SparseCOOSolverCreator.h"
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
#include "../../include/std/EvalNewPointStd_StepSetOptions.h"
#include "../../include/std/EvalNewPointTailoredApproachStd_Step.h"
#include "../../include/std/EvalNewPointTailoredApproachStd_StepSetOptions.h"
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
#include "../../include/std/CheckConvergenceStd_AddedStepSetOptions.h"
#include "../../include/std/CheckSkipBFGSUpdateStd_Step.h"
#include "../../include/std/CheckSkipBFGSUpdateNoPy_Step.h"
#include "../../include/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.h"
#include "../../include/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.h"
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

#include "LinAlgPack/include/PermVecMat.h"

#include "Misc/include/dynamic_cast_verbose.h"

// Stuff for exact reduced hessian
#include "../../include/std/ReducedHessianExactStd_Step.h"
#include "../../include/std/CrossTermExactStd_Step.h"
#include "../../include/std/DampenCrossTermStd_Step.h"
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefChol.h"	// rHL_k

// Correct a bad initial guess
#include "../../include/std/CorrectBadInitGuessStd_AddedStep.h"
#include "../../include/std/CorrectBadInitGuessStd_AddedStepSetOptions.h"

namespace ReducedSpaceSQPPack {

rSQPAlgo_ConfigMamaJama::rSQPAlgo_ConfigMamaJama(
		  ReferenceCountingPack::ref_count_ptr<BasisSystem> basis_sys_ptr
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			Gc_iq_creator_ptr
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			U_iq_creator_ptr
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			HL_iq_creator_ptr
		)
	: basis_sys_ptr_(basis_sys_ptr)
		, Gc_iq_creator_ptr_(Gc_iq_creator_ptr)
		, U_iq_creator_ptr_(U_iq_creator_ptr)
		, HL_iq_creator_ptr_(HL_iq_creator_ptr)
		, qp_solver_type_(QPOPT)
		, linear_solver_type_(MA28)
		, factorization_type_(AUTO_FACT)
		, bigM_(-1.0)
		, correct_bad_init_guess_(false)
		, max_basis_cond_change_frac_(-1.0)
		, warm_start_frac_(-1.0)
		, merit_function_type_(MERIT_FUNC_MOD_L1_INCR)
		, merit_function_penalty_param_update_(MERIT_FUNC_PENALTY_PARAM_WITH_MULT)
		, line_search_method_(LINE_SEARCH_DIRECT)
		, use_line_search_correct_kkt_tol_(-1.0)
		, full_steps_after_k_(-1)
		, exact_reduced_hessian_(false)
		, quasi_newton_(QN_AUTO)
		, hessian_initialization_(INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS)
		, quasi_newton_dampening_(false)
		, num_lbfgs_updates_stored_(7)
		, lbfgs_auto_scaling_(true)
		, max_dof_quasi_newton_dense_(1000)
		, check_results_(false)
{}

rSQPAlgo_ConfigMamaJama::~rSQPAlgo_ConfigMamaJama() {
	// No need to realy do anything!
}

void rSQPAlgo_ConfigMamaJama::warm_start_frac(value_type warm_start_frac) {
	if( warm_start_frac <= 0.0 || warm_start_frac > 1.0 )
		throw std::logic_error( "rSQPAlgo_ConfigMamaJama::warm_start_frac(frac) : "
			"Error, frac must be in the range (0,1]"								);
	warm_start_frac_ = warm_start_frac;
}

void rSQPAlgo_ConfigMamaJama::config_algo_cntr(rSQPAlgoContainer& algo_cntr
	, std::ostream* trase_out)
{

//	TRACE0( "\n*** rSQPAlgo_ConfigMamaJama::config_algo_cntr(algo_cntr) called\n" );
//	TRACE0( "\n*** Before algo is configured\n" );
//	print_state();

	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	
	// ////////////////////////////////////////////////////////////
	// A. ???

	// /////////////////////////////////////////////////////////////////////////
	// B. Create an algo object, give to algo_cntr, then give algo_cntr to algo

#ifdef RELEASE_TRACE
	std::cout << "\n*** Create the algo object ...\n";
#endif

	namespace rcp = ReferenceCountingPack;
	typedef rcp::ref_count_ptr<rSQPAlgo>	algo_ptr_t;
	
	algo_ptr_t algo(new rSQPAlgo);
	assert(algo.get());

#ifdef RELEASE_TRACE
	std::cout << "\n*** algo_cntr.set_algo(algo.get()) ...\n";
#endif

	algo_cntr.set_algo( rcp::rcp_implicit_cast<rSQPAlgoInterface>(algo) );

#ifdef RELEASE_TRACE
	std::cout << "\n*** algo->set_algo_cntr(&algo_cntr) ...\n";
#endif

	algo->set_algo_cntr(&algo_cntr);

	// /////////////////////////////////////////////
	// C. Configure algo

	// /////////////////////////////////////////////////////
	// C.0 Set the nlp and track objects

#ifdef RELEASE_TRACE
	std::cout << "\n*** Set the NLP and track objects ...\n";
#endif

	algo->set_nlp( algo_cntr.get_nlp().get() );
	algo->set_track( rcp::rcp_implicit_cast<AlgorithmTrack>(algo_cntr.get_track()) );

	// 7/28/00: Determine if this is a standard NLPReduced nlp or a tailored approach nlp.

	bool tailored_approach
		= (NULL != dynamic_cast<NLPrSQPTailoredApproach*>(algo->get_nlp()));
	if( tailored_approach ) {
		// Change the options for the tailored approach. 

		if(trase_out) {
			*trase_out
				<< "\n***********************************************************************************\n"
				<< "This is a tailored approach NLP and the following options are used or not allowed:\n"
				<< "merit_function_type = L1;\n"
				<< "merit_function_penalty_param_update = MULT_FREE;\n"
				<< "qp_solver != QPSCPD; *** set to QPKWIK if QPSCPD is used\n"
				<< "fact_type_ = DIRECT_FACT;\n"
				;
		}

		merit_function_type_					= MERIT_FUNC_L1;
		merit_function_penalty_param_update_	= MERIT_FUNC_PENALTY_PARAM_MULT_FREE;
		if( qp_solver_type_ == QPSCPD )
			qp_solver_type_ = QPKWIK;
		factorization_type_ = DIRECT_FACT;

	}
	else {
		if( NULL != dynamic_cast<NLPrSQPTailoredApproach*>(algo->get_nlp()) ) {
			std::ostringstream omsg;
			omsg
				<< "rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : "
				<< "Error, type nlp object with the concrete type "
				<< typeid(algo->nlp()).name()
				<< " does not support the NLPReduced interface.";
			if(trase_out)
				*trase_out << omsg.str() << std::endl;
			throw std::logic_error( omsg.str() );
		}
	}

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
			if( ! algo->nlp().has_bounds()  || qp_solver_type_ == VE09
					|| qp_solver_type_ == QPSCPD )
			{
				fact_type_ = ADJOINT_FACT;
			}
			else {
				// If the total number of bounds is less than the number
				// of degrees of freedom then the adjoint factorization
				// will be faster.
				NLP &nlp = algo->nlp();
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

		matrix_iqa_creator_ptr_t
			matrix_iqa_creators[matrix_creator_t::num_MatrixWithOp];
		matrix_iqa_creators[matrix_creator_t::Q_HL]
			= HL_iq_creator_ptr_.get()
				? HL_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>();
		matrix_iqa_creators[matrix_creator_t::Q_Gc]
			= Gc_iq_creator_ptr_.get()
				? Gc_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<COOMatrixWithPartitionedViewSubclass>();
		matrix_iqa_creators[matrix_creator_t::Q_Y]
			= new IterQuantMatrixWithOpCreatorContinuous<IdentZeroVertConcatMatrixSubclass>();
		matrix_iqa_creators[matrix_creator_t::Q_Z]
			= NULL;		// Determine later
		matrix_iqa_creators[matrix_creator_t::Q_U]
			= U_iq_creator_ptr_.get()
				? U_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<COOMatrixPartitionViewSubclass>();
		matrix_iqa_creators[matrix_creator_t::Q_V]
			= new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>();
		matrix_iqa_creators[matrix_creator_t::Q_rHL]
			= NULL;		// Determine later

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
			n = algo->nlp().n(),
			m = algo->nlp().m();
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

		if( ! algo->nlp().has_bounds() ) {
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

		storage_num[state_t::Q_kkt_err]			= 2;
		storage_num[state_t::Q_rGL]				= 2;
		storage_num[state_t::Q_lambda]			= 2;

		storage_num[state_t::Q_Ypy]				= 2;
		storage_num[state_t::Q_Zpz]				= 2;

		storage_num[state_t::Q_d]				= 2;

		storage_num[state_t::Q_rGf]				= 2;

		storage_num[state_t::Q_alpha]			= 2;
		storage_num[state_t::Q_mu]				= 2;
		storage_num[state_t::Q_phi]				= 2;
		storage_num[state_t::Q_nu]				= 2;

		// create a new state object
		algo->set_state( new state_t( new matrix_creator_aggr_t(matrix_iqa_creators), storage_num ) );

		// Now set the number of LBFGS update vectors
		if(use_limited_memory) {
			// Here we will set the number of update vectors to store by setting
			// it for the k-1 iteration.  That way the k iteration will still not
			// be updated and therefore the regular initializations will still
			// be performed.  This is a little bit of a hack but it should work.
			SymLBFGSMatrixSubclass *_rHL
				= dynamic_cast<SymLBFGSMatrixSubclass*>(&algo->rsqp_state().rHL().set_k(-1));
			if(_rHL) {
				_rHL->set_num_updates_stored( num_lbfgs_updates_stored_ );
				_rHL->auto_rescaling( lbfgs_auto_scaling_ );
				algo->rsqp_state().rHL().set_not_updated(-1);
			}
		}
	}

	// /////////////////////////////////////////////////////
	// C.2. Create and set the decomposition system object

#ifdef RELEASE_TRACE
	std::cout << "\n*** Create and set the decomposition system object ...\n";
#endif

	DecompositionSystemVarReductStd*
		decomp_sys = NULL;

	if( !tailored_approach ) {
	
		// Only need a decomposition system object if we are not using the tailored approach

		if( !basis_sys_ptr_.get() ) {
			SparseCOOSolverCreator
				*sparse_solver_creator = NULL;
			if( linear_solver_type_ == MA28 ) {
				sparse_solver_creator
					= new SparseSolverPack::MA28SparseCOOSolverCreator(
							new SparseSolverPack::MA28SparseCOOSolverSetOptions
							, const_cast<OptionsFromStreamPack::OptionsFromStream*>(options_)
						);
			}
			else if ( linear_solver_type_ == MA48 ) {
				sparse_solver_creator
					= new SparseSolverPack::MA48SparseCOOSolverCreator(
							new SparseSolverPack::MA48SparseCOOSolverSetOptions
							, const_cast<OptionsFromStreamPack::OptionsFromStream*>(options_)
						);
			}
			else {
				assert(0);
			}

			// Note, the switch based code for the above if statements would not
			// compile under MipsPro 7.3.1.1m.

			basis_sys_ptr_ = new COOBasisSystem(sparse_solver_creator, true);
		}

		DecompositionSystemVarReductImpNode* decomp_sys_aggr = 0;
		
		if ( fact_type_ == DIRECT_FACT ) {
			decomp_sys_aggr = new DecompositionSystemCoordinateDirect(
									basis_sys_ptr_.release(), true);
		}
		else if ( fact_type_ == ADJOINT_FACT ) {
			decomp_sys_aggr = new DecompositionSystemCoordinateAdjoint(
									basis_sys_ptr_.release(), true);
		}
		else {
			assert(0);
		}

		// Note, the switch based code for the above if statements would not
		// compile under MipsPro 7.3.1.1m.

		DecompositionSystemVarReductStd*
			decomp_sys = new DecompositionSystemVarReductStd( false, algo.get(), decomp_sys_aggr );

		algo->rsqp_state().set_decomp_sys(decomp_sys);

	}


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
		if( tailored_approach ) {
			typedef rcp::ref_count_ptr<NLPrSQPTailoredApproachTester>
				nonconst_deriv_tester_ptr_t;
			typedef EvalNewPointTailoredApproachStd_Step::deriv_tester_ptr_t
				deriv_tester_ptr_t;
			
			nonconst_deriv_tester_ptr_t
				deriv_tester = new NLPrSQPTailoredApproachTester();

			{
				NLPrSQPTailoredApproachTesterSetOptions options_setter(deriv_tester.get());
				options_setter.set_options(*options_);
			}		

			EvalNewPointTailoredApproachStd_Step
				*eval_new_point_step = new EvalNewPointTailoredApproachStd_Step(deriv_tester);

			{
				EvalNewPointTailoredApproachStd_StepSetOptions options_setter(eval_new_point_step);
				options_setter.set_options(*options_);
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
				options_setter.set_options(*options_);
			}		

			EvalNewPointStd_Step
				*eval_new_point_step = new EvalNewPointStd_Step(deriv_tester);

			{
				EvalNewPointStd_StepSetOptions options_setter(eval_new_point_step);
				options_setter.set_options(*options_);
			}		

			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );
		}

		// (2) ReducedGradient
		algo->insert_step( ++step_num, ReducedGradient_name, new ReducedGradientStd_Step );

		// (3) Calculate Reduced Gradient of the Lagrangian
		algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, new CalcReducedGradLagrangianStd_AddedStep );

		// (4)	Calculate the Lagrange multipliers for the independent constraints.
		// 		These are computed here just in case the algorithm converges and we need to
		// 		report these multipliers to the NLP.
		if( !tailored_approach ) {
			algo->insert_step( ++step_num, CalcLambdaIndep_name, new  CalcLambdaIndepStd_AddedStep );
		}

		// (5) Check for convergence
		{
			CheckConvergenceStd_AddedStep
				*check_convergence_step = new CheckConvergenceStd_AddedStep;

			CheckConvergenceStd_AddedStepSetOptions
				opt_setter( check_convergence_step );
			opt_setter.set_options( *options_ );

			algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
		}


		// (6) ReducedHessian
		ref_count_ptr<ReducedHessianBFGSStd_Step>
			bfgs_updater = new ReducedHessianBFGSStd_Step(quasi_newton_dampening_);
		algo->insert_step( ++step_num, ReducedHessian_name
			, rcp::rcp_implicit_cast<AlgorithmStep>(bfgs_updater) );

		// (6.-1) CheckSkipBFGSUpdate
		algo->insert_assoc_step(
			  step_num
			, GeneralIterationPack::PRE_STEP
			, 1
			, CheckSkipBFGSUpdate_name
			, new CheckSkipBFGSUpdateStd_Step
		  );

		// (7) DepDirec
		if( !tailored_approach ) {
			algo->insert_step( ++step_num, DepDirec_name, new  DepDirecStd_Step );
		}

		// (8) IndepDirec
		algo->insert_step( ++step_num, IndepDirec_name, new  IndepDirecWithoutBounds_Step );

		// (9) CalcDFromYPYZPZ

		algo->insert_step( ++step_num, "CalcDFromYPYZPZ", new CalcDFromYPYZPZ_Step );

		// (10) LineSearch

		if( line_search_method_ == LINE_SEARCH_NONE ) {

			algo->insert_step(
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
					switch( algo_cntr.journal_output_level() ) {
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

			algo->insert_step(
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
				case MERIT_FUNC_L1: {
					switch(merit_function_penalty_param_update_) {
						case MERIT_FUNC_PENALTY_PARAM_WITH_MULT:
							param_update_step
								= new  MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
											rcp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
							break;
						case MERIT_FUNC_PENALTY_PARAM_MULT_FREE:
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
			ppu_options_setter.set_options( *options_ );

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

				algo->insert_assoc_step(	  step_num
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

		if( algo_cntr.nlp().has_bounds() ) {
			
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
					if(algo_cntr.journal_output_level() == PRINT_ITERATION_QUANTITIES)
						mapped_qp_file_ = mapped_qp_file_ptr_t( new std::ofstream("mapped_qp.txt") );
					_qp_solver->qp_mapping_output( mapped_qp_file_.get() );
					qp_solver = qp_solver_ptr_t( new ReducedQPSolverQPOPTSOLStd(_qp_solver,algo.get()) );
					break;
				}
				case QPKWIK:
				{
					ReducedQPSolverQPKWIKNEW*
						_qp_solver = new ReducedQPSolverQPKWIKNEW;
					if(algo_cntr.journal_output_level() == PRINT_ITERATION_QUANTITIES)
						_qp_solver->create_qpkwiknew_file(true);
					ReducedQPSolverQPKWIKNEWStd*
						__qp_solver = new ReducedQPSolverQPKWIKNEWStd(_qp_solver,algo.get());
					if( warm_start_frac_ > 0.0 )
						__qp_solver->warm_start_frac(warm_start_frac_);
					qp_solver = qp_solver_ptr_t( __qp_solver );
					break;
				}
			}

			Algorithm::poss_type poss;

			poss = algo->get_step_poss(IndepDirec_name);
			algo->replace_step(
						poss
						,new  IndepDirecExact_Step( 
							new ReducedQPSolverCheckOptimality(
								qp_solver, algo.get(), check_results_ )
						  )
				    );

			Algorithm::step_ptr_t calc_rgrad_lagr, calc_lambda, check_conv;

			// Remove and save CalcReducedGradLagr..., CalcLambdaIndep... and CheckConv...

			poss			= algo->get_step_poss(CalcReducedGradLagrangian_name);
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
				*qp_solver = new QPMixedFullReducedQPSCPDSolver(algo.get());
			qp_solver->set_pick_violated_policy(
					ConstrainedOptimizationPack::QPSCPDPack::ConstraintsVarBoundsRelaxed::MOST_VIOLATED );
			if( bigM_ > 0.0 )
				qp_solver->set_M(bigM_);
			if( warm_start_frac_ > 0.0 )
				qp_solver->set_warm_start_frac(warm_start_frac_);
			using ConstrainedOptimizationPack::QPSCPD;	
			QPSCPD::EOutputLevel	qp_olevel;
			switch( algo_cntr.journal_output_level() ) {
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
					 qp_solver, algo.get(), check_results_ );

			Algorithm::poss_type poss;
			poss = algo->get_step_poss(IndepDirec_name);

			algo->remove_step(poss);

			algo->insert_step(
				poss
				, SearchDirec_name
				, new SearchDirecMixedFullReduced_Step(
					qp_opt_checker
				  )
			  );

			algo->insert_step(
				poss + 1
				, "CalcZpzFromDYPY"
				, new CalcZpzFromDYPY_Step
			  );

			poss = algo->get_step_poss(CalcLambdaIndep_name);
			algo->remove_step(poss);

		}
	}

	// 8/30/99: Add the iteration quantities for the QPSolverStats and 
	// ActSetStats and the added step that will calculate the
	// active set updates.

	{
		// Add active set and QP statistics to state object
		algo->state().set_iter_quant( act_set_stats_name
			, new IterQuantityAccessContinuous<ActSetStats>( 1, act_set_stats_name ) );
		algo->state().set_iter_quant( qp_solver_stats_name
			, new IterQuantityAccessContinuous<QPSolverStats>( 1, qp_solver_stats_name ) );

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

		// Output of KKT conditions of reduced QP subproblem
		algo->insert_assoc_step( poss, GeneralIterationPack::POST_STEP, nas+3
			, "PrintReducedQPError", new PrintReducedQPError_AddedStep( print_qp_error_ ) );

	}

	// 10/15/99: Add basis checking step
	if( !tailored_approach ) {
		CheckBasisFromPy_Step
			*basis_check_step = new CheckBasisFromPy_Step(
				new NewBasisSelectionStd_Strategy(decomp_sys) );
		if( max_basis_cond_change_frac_ > 0.0 )
			basis_check_step->max_basis_cond_change_frac( max_basis_cond_change_frac_ );
		Algorithm::poss_type poss;
		assert(poss = algo->get_step_poss( DepDirec_name ) );
		algo->insert_step(
			  poss+1
			, "CheckBasis"
			, basis_check_step
		  );
	}

	// 10/19/99: Add quasi-newton stats
	algo->state().set_iter_quant( quasi_newton_stats_name
		, new IterQuantityAccessContinuous<QuasiNewtonStats>( 1, quasi_newton_stats_name ) );

	// 11/12/99: Add the option of using full steps after a specified number of iterations
	// if this option has been set and if we are using a linesearch method.
	if( full_steps_after_k_ > 0 && line_search_method_ != LINE_SEARCH_NONE ) {
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
				, full_steps_after_k_
			  )
		  );

	}

	// 12/3/99: Adding finite difference initializaiton of the reduced hessian
	if( hessian_initialization_ != INIT_HESS_IDENTITY ) {
		Algorithm::poss_type poss;
		assert(poss = algo->get_step_poss( ReducedHessian_name ) );
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
		algo->insert_assoc_step( poss, GeneralIterationPack::PRE_STEP, 1
			, "InitFiniteDiffReducedHessian"
			, init_red_hess_step  );
	}

	// 6/13/00: Adding steps for the computation of the exact reduced Hessian
	if( exact_reduced_hessian_ == true ) {

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
		DampenCrossTermStd_Step
			*zeta_step = new DampenCrossTermStd_Step;
		// ToDo: set options from stream
		algo->insert_assoc_step( poss+1, GeneralIterationPack::POST_STEP, 1
			, "DampenReducedQPCrossTerm"
			, zeta_step  );

		// Change the type of the iteration quantity for rHL
		typedef GeneralIterationPack::IterQuantityAccessContinuous<
					ConstrainedOptimizationPack::MatrixSymPosDefChol >
				iq_rHL_concrete_t;
		typedef GeneralIterationPack::IterQuantityAccessDerivedToBase< MatrixWithOp
					, ConstrainedOptimizationPack::MatrixSymPosDefChol >
				iq_rHL_t;
		algo->state().get_iter_quant(algo->state().get_iter_quant_id(rHL_name))
			= AlgorithmState::IQ_ptr( new iq_rHL_t( new iq_rHL_concrete_t(1,rHL_name) ) );

		// That's it man!
	}

	// 7/3/00: Adding some steps to correct a bad initial guess.
	if( correct_bad_init_guess_ ) {
		typedef rcp::ref_count_ptr<CorrectBadInitGuessStd_AddedStep>
			corr_xinit_ptr_t;
		corr_xinit_ptr_t
			corr_xinit_ptr = new CorrectBadInitGuessStd_AddedStep;
		CorrectBadInitGuessStd_AddedStepSetOptions
			options_setter(corr_xinit_ptr.get());
		options_setter.set_options(*options_);

		Algorithm::poss_type poss;

		// Add it after the computation of the reduced gradient of the Lagrangian
		assert(poss = algo->get_step_poss( CalcReducedGradLagrangian_name ) );
		algo->insert_step(
			  poss+1
			, "CorrectBadInitGuessFirst"
			, rcp::rcp_implicit_cast<AlgorithmStep>( corr_xinit_ptr )
		  );

		// Add it after the line search
		assert(poss = algo->get_step_poss( LineSearch_name ) );
		algo->insert_step(
			  poss+1
			, "CorrectBadInitGuessSecond"
			, rcp::rcp_implicit_cast<AlgorithmStep>( corr_xinit_ptr )
		  );
	}

// Desprite debugging stuff
//	assert(algo.get());
//	TRACE0( "\n*** After algo is configured\n" );
//	print_state();

}

void rSQPAlgo_ConfigMamaJama::init_algo(rSQPAlgoInterface& _algo)
{

//	TRACE0( "\n*** rSQPAlgo_ConfigMamaJama::init_algo(_algo) called\n" );
//	TRACE0( "\n*** Before algo is initialized\n" );
//	print_state();

	using DynamicCastHelperPack::dyn_cast;

	namespace rcp = ReferenceCountingPack;

	rSQPAlgo	&algo	= dyn_cast<rSQPAlgo>(_algo);
	rSQPState	&state	= algo.rsqp_state();
	NLP			&nlp = algo.nlp();

	algo.max_iter( algo.algo_cntr().max_iter() );
	algo.rsqp_state().check_results( check_results_ );

	// (1) Initialize the nlp
	nlp.initialize();

	// Advance everything two iterations to wipe out any memory to k-1
	state.next_iteration();
	state.next_iteration();

	// (2) Set the initial x to the initial guess.
	state.x().set_k(0).v() = nlp.xinit();

	// (3) Get the initial basis from the nlp.

	const size_type
		n = nlp.n(),
		m = nlp.m(),
		r = nlp.r();
	state.var_dep()		= Range1D(1, r);
	state.var_indep()	= Range1D(r + 1, n);
	state.con_indep()	= Range1D(1, r);
	state.con_dep()		= (r < m) ? Range1D(r + 1, m) : Range1D(m + 1, m + 1);

	if( NLPReduced	*red_nlp = dynamic_cast<NLPReduced*>(algo.get_nlp()) ) {
		// The nlp will have a basis
		// set whether or not nlp.nlp_selects_basis() returns true (See NLPReduced).
		size_type rank;
		red_nlp->get_basis( &state.var_perm_new(), &state.con_perm_new(), &rank );
	}
	else if ( dynamic_cast<NLPrSQPTailoredApproach*>(algo.get_nlp()) ){
		// Just set to original order (identity).
		using LinAlgPack::identity_perm;
		state.var_perm_new().resize(n);
		identity_perm(&state.var_perm_new());
		state.con_perm_new().resize(m);
		identity_perm(&state.con_perm_new());
	}

	state.initialize_fast_access();

	// Reset the iteration count to zero
	state.k(0);

	// Get organized output of vectors and matrices even if setw is not used.
	algo.track().journal_out()
		<< std::setprecision(algo.algo_cntr().print_digits())
		<< std::scientific;

	// set the first step
	algo.do_step_first(1);

	// The rest of the algorithm will initialize itself
}

}	// end namespace ReducedSpaceSQPPack 