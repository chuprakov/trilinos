// /////////////////////////////////////////////////////////////////////////
// rSQPAlgo_ConfigMamaJama.cpp

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
#include "ReducedSpaceSQPPack/include/rSQPStateContinuousStorage.h"
#include "ReducedSpaceSQPPack/include/rSQPStateContinuousStorageMatrixWithOpCreatorAggr.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/COOMatrixWithPartitionedViewSubclass.h"			// HL, Gc, Hcj
#include "ConstrainedOptimizationPack/include/IdentZeroVertConcatMatrixSubclass.h"	// Y
#include "ConstrainedOptimizationPack/include/DenseIdentVertConcatMatrixSubclass.h"	// Z
#include "ConstrainedOptimizationPack/include/ZAdjointFactMatrixSubclass.h"			// .
#include "SparseLinAlgPack/include/COOMatrixPartitionViewSubclass.h"				// U
#include "SparseLinAlgPack/include/GenMatrixSubclass.h"								// V
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefCholFactor.h"          // rHL 
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefInvCholFactor.h"		// .
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefLBFGS.h"				// .
#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasicInitDiagonal.h"// | rHL (super basics)
#include "SparseLinAlgPack/include/MatrixSymDiagonalStd.h"                          // |

#include "ConstrainedOptimizationPack/include/VariableBoundsTesterSetOptions.h"

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
#include "ReducedSpaceSQPPack/include/std/DecompositionSystemVarReductStd.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateDirect.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateAdjoint.h"
#include "ConstrainedOptimizationPack/include/DirectLineSearchArmQuad_StrategySetOptions.h"

#include "SparseSolverPack/test/BasisSystemTesterSetOptions.h"

#include "ConstrainedOptimizationPack/include/QPSolverRelaxedTester.h"
#include "ConstrainedOptimizationPack/include/QPSolverRelaxedTesterSetOptions.h"
#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPSchur.h"
#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPSchurSetOptions.h"
#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianFull.h"
#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianSuperBasic.h"
#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPKWIK.h"

#include "ReducedSpaceSQPPack/include/std/ReducedQPSolverCheckOptimality.h"
#include "ReducedSpaceSQPPack/include/std/ReducedQPSolverQPOPTSOLStd.h"
#include "ReducedSpaceSQPPack/include/std/ReducedQPSolverQPSOL.h"
#include "ReducedSpaceSQPPack/include/std/ReducedQPSolverQPOPT.h"

#include "ReducedSpaceSQPPack/include/std/rSQPAlgorithmStepNames.h"

#include "ReducedSpaceSQPPack/include/std/DecompositionSystemVarReductStd.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproach_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachCoordinate_Step.h"
#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachOrthogonal_Step.h"
#include "ReducedSpaceSQPPack/include/std/ReducedGradientStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/InitFinDiffReducedHessian_Step.h"
#include "ReducedSpaceSQPPack/include/std/InitFinDiffReducedHessian_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSFull_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSProjected_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateLPBFGS_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.h"
#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_StrategySetOptions.h"
#include "ReducedSpaceSQPPack/include/std/DepDirecStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/CheckBasisFromCPy_Step.h"
#include "ReducedSpaceSQPPack/include/std/CheckBasisFromPy_Step.h"
#include "ReducedSpaceSQPPack/include/std/IndepDirecWithoutBounds_Step.h"
#include "ReducedSpaceSQPPack/include/std/SetDBoundsStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/IndepDirecWithBoundsStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/IndepDirecWithBoundsStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/IndepDirecExact_Step.h"
#include "ReducedSpaceSQPPack/include/std/CalcDFromYPYZPZ_Step.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchFailureNewBasisSelection_Step.h"
#include "ReducedSpaceSQPPack/include/std/NewBasisSelectionStd_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchFullStep_Step.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchDirect_Step.h"
#include "ReducedSpaceSQPPack/include/std/LineSearch2ndOrderCorrect_Step.h"
#include "ReducedSpaceSQPPack/include/std/LineSearch2ndOrderCorrect_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchWatchDog_Step.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchWatchDog_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/LineSearchFullStepAfterKIter_Step.h"
#include "ReducedSpaceSQPPack/include/std/CalcLambdaIndepStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CalcReducedGradLagrangianStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CheckConvergenceStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CheckConvergenceStd_AddedStepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/CheckSkipBFGSUpdateStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/CheckSkipBFGSUpdateStd_StepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamUpdateWithMult_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_ModifiedL1LargerSteps_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.h"
#include "ReducedSpaceSQPPack/include/std/ActSetStats_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/NumFixedDepIndep_AddedStep.h"

#include "ReducedSpaceSQPPack/include/std/act_set_stats.h"
#include "ReducedSpaceSQPPack/include/std/qp_solver_stats.h"
#include "ReducedSpaceSQPPack/include/std/quasi_newton_stats.h"

#include "ConstrainedOptimizationPack/include/DecompositionSystemCoordinateDirect.h"
#include "SparseLinAlgPack/include/sparse_bounds.h"

#include "LinAlgPack/include/PermVecMat.h"

#include "Misc/include/dynamic_cast_verbose.h"
#include "Misc/include/ReleaseResource_ref_count_ptr.h"

// Stuff to readin options
#include "Misc/include/StringToIntMap.h"
#include "Misc/include/StringToBool.h"

// Stuff for exact reduced hessian
#include "ReducedSpaceSQPPack/include/std/ReducedHessianExactStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/CrossTermExactStd_Step.h"
#include "ReducedSpaceSQPPack/include/std/DampenCrossTermStd_Step.h"

// Correct a bad initial guess
#include "ReducedSpaceSQPPack/include/std/CorrectBadInitGuessStd_AddedStep.h"
#include "ReducedSpaceSQPPack/include/std/CorrectBadInitGuessStd_AddedStepSetOptions.h"

namespace {
	const int INF_BASIS_COND_CHANGE_FRAC         = 1e+20;
	const int DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE = 200;
}

namespace ReducedSpaceSQPPack {

//
// Here is where we defint the default values for the algorithm.  These
// should agree with what are in the rSQPpp.opt.rsqp_mama_jama_solve file.
//
rSQPAlgo_ConfigMamaJama::SOptionValues::SOptionValues()
	:
		  direct_linear_solver_type_(LA_AUTO)
		, null_space_matrix_type_(NULL_SPACE_MATRIX_AUTO)
		, range_space_matrix_type_(RANGE_SPACE_MATRIX_AUTO)
		, max_basis_cond_change_frac_(-1.0)
		, exact_reduced_hessian_(false)
		, quasi_newton_(QN_AUTO)
		, max_dof_quasi_newton_dense_(-1)
		, num_lbfgs_updates_stored_(-1)
		, lbfgs_auto_scaling_(true)
		, hessian_initialization_(INIT_HESS_AUTO)
		, qp_solver_type_(QP_AUTO)
		, line_search_method_(LINE_SEARCH_AUTO)
		, merit_function_type_(MERIT_FUNC_AUTO)
		, l1_penalty_param_update_(L1_PENALTY_PARAM_AUTO)
		, full_steps_after_k_(-1)
		, print_qp_error_(false)
		, correct_bad_init_guess_(false)
{}

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
		, options_(NULL)
{}

rSQPAlgo_ConfigMamaJama::~rSQPAlgo_ConfigMamaJama() {
	// No need to really do anything!
}

void rSQPAlgo_ConfigMamaJama::set_options( const OptionsFromStreamPack::OptionsFromStream* options )
{
	options_ = options;
}

void rSQPAlgo_ConfigMamaJama::config_algo_cntr(rSQPAlgoContainer& algo_cntr
	, std::ostream* trase_out)
{
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;

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
	algo_ptr_t algo(new rSQPAlgo);
	assert(algo.get());
	algo_cntr.set_algo( rcp::rcp_implicit_cast<rSQPAlgoInterface>(algo) );
	algo->set_algo_cntr(&algo_cntr);

	// /////////////////////////////////////////////
	// C. Configure algo

	// /////////////////////////////////////////////////////
	// C.0 Set the nlp and track objects

	if(trase_out)
		*trase_out << "\nSetting the NLP and track objects ...\n";

	algo->set_nlp( algo_cntr.get_nlp().get() );
	algo->set_track( rcp::rcp_implicit_cast<AlgorithmTrack>(algo_cntr.get_track()) );

	// ////////////////////////////////////////////////
	// Determine what the options are:

	// Readin the options
	if(options_) {
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
	using SparseLinAlgPack::num_bounds;
	const size_type
		n = algo->nlp().n(),
		r = algo->nlp().r(),
		dof = nlp.n() - nlp.r(),
		nb = nlp.has_bounds() ? num_bounds( nlp.xl(), nlp.xu() ) : 0;

	// 7/28/00: Determine if this is a standard NLPReduced nlp or a tailored approach nlp.
	bool tailored_approach
		= (NULL != dynamic_cast<NLPrSQPTailoredApproach*>(algo->get_nlp()));
	if( tailored_approach ) {
		// Change the options for the tailored approach. 
		if(trase_out) {
			*trase_out
				<< "\nThis is a tailored approach NLP and the following options are set:\n"
				<< "merit_function_type         = L1;\n"
				<< "l1_penalty_parameter_update = MULT_FREE;\n"
				<< "null_space_matrix           = EXPLICIT;\n"
				;
		}
		cov_.merit_function_type_		= MERIT_FUNC_L1;
		cov_.l1_penalty_param_update_	= L1_PENALTY_PARAM_MULT_FREE;
		cov_.null_space_matrix_type_    = NULL_SPACE_MATRIX_EXPLICIT;
	}
	else {
		// If it is not a tailored approach NLP then it better
		// support the NLPReduced interface!
		if( NULL == dynamic_cast<NLPReduced*>(algo->get_nlp()) ) {
			std::ostringstream omsg;
			omsg
				<< "rSQPAlgo_ConfigMamaJama::config_algo_cntr(...) : "
				<< "Error, type nlp object with the concrete type "
				<< typeid(algo->nlp()).name()
				<< " does not support the NLPReduced interface.";
			if(trase_out)
				*trase_out << std::endl << omsg.str() << std::endl;
			throw std::logic_error( omsg.str() );
		}
	}

	// Determine whether to use direct or adjoint factorization
	switch(uov_.null_space_matrix_type_) {
		case NULL_SPACE_MATRIX_EXPLICIT:
			cov_.null_space_matrix_type_ = NULL_SPACE_MATRIX_EXPLICIT;
			break;
		case NULL_SPACE_MATRIX_IMPLICIT:
			cov_.null_space_matrix_type_ = NULL_SPACE_MATRIX_IMPLICIT;
			break;
		case NULL_SPACE_MATRIX_AUTO:
		{
			if( !nlp.has_bounds() || uov_.qp_solver_type_ == QP_QPSCHUR )
			{
				if(trase_out) {
					*trase_out << "\nnull_sapce_matrix == AUTO:";
					if(!nlp.has_bounds())
						*trase_out << "\nnlp.has_bounds() == " << algo->nlp().has_bounds() << ":";
					else if( uov_.qp_solver_type_ == QP_QPSCHUR )
						*trase_out << "\nqp_solver == QPSCHUR:";
					else
						assert(0);
					*trase_out << "\nsetting null_space_matrix = IMPLICIT ...\n";
				}
				cov_.null_space_matrix_type_ = NULL_SPACE_MATRIX_IMPLICIT;
			}
			else {
				if(trase_out)
					*trase_out << "\nnull_sapce_matrix == AUTO:\n"
									"Checking the total number of variable bounds ...\n";
				// If the total number of bounds is less than the number
				// of degrees of freedom then the adjoint factorization
				// will be faster.
				if(trase_out)
					*trase_out << "#bounds = " << nb << ", n-r = " << dof << std::endl;
				if( nb < dof ) { 
					if(trase_out)
						*trase_out
							<< "There are fewer bounds than degrees of freedom:\n"
							"setting null_space_matrix = IMPLICIT ...\n";
					cov_.null_space_matrix_type_ = NULL_SPACE_MATRIX_IMPLICIT;
				}
				else {
					if(trase_out)
						*trase_out
							<< "There are more bounds than degrees of freedom:\n"
							"setting null_space_matrix = EXPLICIT ...\n";
					cov_.null_space_matrix_type_ = NULL_SPACE_MATRIX_EXPLICIT;
				}
			}
			break;
		}
	}

	// ToDo: Implement the orthogonal decompositon for general NLPs
	if( !tailored_approach
		&& uov_.range_space_matrix_type_ != RANGE_SPACE_MATRIX_COORDINATE )
	{
		if(trase_out)
			*trase_out << "\nrange_space_matrix != COORDINATE:\n"
							"Sorry, the orthogonal decomposition is not "
							"supported for general NLPs yet!\n"
							"setting range_space_matrix = COORDINATE ...\n";
		cov_.range_space_matrix_type_ = RANGE_SPACE_MATRIX_COORDINATE;
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
					<< "\nquasi_newton == AUTO:\n"
					<< "\nnlp.has_bounds() == " << nlp.has_bounds() << ":\n";
			if( n - r > cov_.max_dof_quasi_newton_dense_ ) {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " > max_dof_quasi_newton_dense = "
						<< cov_.max_dof_quasi_newton_dense_ <<  ":\n";
				if( !algo->nlp().has_bounds() ) {
					if(trase_out)
						*trase_out
							<< "setting quasi_newton == LPBFGS\n";
					cov_.quasi_newton_ = QN_LPBFGS;
				}
				else {
					if(trase_out)
						*trase_out
							<< "setting quasi_newton == LBFGS\n";
					cov_.quasi_newton_ = QN_LPBFGS;
				}
			}
			else {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " <= max_dof_quasi_newton_dense = "
						<< cov_.max_dof_quasi_newton_dense_ << ":\n";
				if( !algo->nlp().has_bounds() ) {
					if(trase_out)
						*trase_out
							<< "setting quasi_newton == BFGS\n";
					cov_.quasi_newton_ = QN_BFGS;
				}
				else {
					if(trase_out)
						*trase_out
							<< "setting quasi_newton == PBFGS\n";
					cov_.quasi_newton_ = QN_PBFGS;
				}
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

	// Set the default options that where not already set yet
	set_default_options(uov_,&cov_,trase_out);

	// Adjust quasi-newton for QPKWIK
	switch( cov_.quasi_newton_ ) {
		case QN_BFGS:
			break;
		case QN_PBFGS:
		case QN_LBFGS:
		case QN_LPBFGS:
			if( cov_.qp_solver_type_ == QP_QPKWIK ) {
				if(trase_out)
					*trase_out
						<< "\nqp_solver == QPKWIK and quasi_newton != BFGS:\n"
						<< "QPKWIK can not use any other option than BFGS, setting quasi_newton = BFGS\n";
				cov_.quasi_newton_ = QN_BFGS;
			}
			break;
	    default:
			assert(0); // Invalid option!
	}

	// /////////////////////////////////////////////////////
	// C.1. Create and set the state object

	if(trase_out)
		*trase_out << "\nCreating and setting the state object ...\n";

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
		// HL
		matrix_iqa_creators[matrix_creator_t::Q_HL]
			= HL_iq_creator_ptr_.get()
				? HL_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>();
		// Gc
		matrix_iqa_creators[matrix_creator_t::Q_Gc]
			= Gc_iq_creator_ptr_.get()
				? Gc_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<COOMatrixWithPartitionedViewSubclass>();
		// Y
		matrix_iqa_creators[matrix_creator_t::Q_Y]
			= new IterQuantMatrixWithOpCreatorContinuous<IdentZeroVertConcatMatrixSubclass>();
		// Z
		switch(cov_.null_space_matrix_type_) {
			case NULL_SPACE_MATRIX_EXPLICIT:
				matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_Z]
					= new IterQuantMatrixWithOpCreatorContinuous<DenseIdentVertConcatMatrixSubclass>();
				break;
			case NULL_SPACE_MATRIX_IMPLICIT:
				matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_Z]
					= new IterQuantMatrixWithOpCreatorContinuous<ZAdjointFactMatrixSubclass>();
				break;
			default:
				assert(0);
		}
		// U
		matrix_iqa_creators[matrix_creator_t::Q_U]
			= U_iq_creator_ptr_.get()
				? U_iq_creator_ptr_
				: new IterQuantMatrixWithOpCreatorContinuous<COOMatrixPartitionViewSubclass>();
		// V
		matrix_iqa_creators[matrix_creator_t::Q_V]
			= new IterQuantMatrixWithOpCreatorContinuous<GenMatrixSubclass>();
		// rHL
		if( ! algo->nlp().has_bounds() ) {
			// Here we just need to solve for linear systems with rHL
			// and we can use a limited memory method or update the dense
			// factors directly.
			matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
				= cov_.quasi_newton_==QN_LBFGS
				? static_cast<IterQuantMatrixWithOpCreator*>(
					new IterQuantMatrixWithOpCreatorContinuous<
					MatrixSymPosDefLBFGS>() )
				: static_cast<IterQuantMatrixWithOpCreator*>( 
					new IterQuantMatrixWithOpCreatorContinuous<
					MatrixSymPosDefInvCholFactor>() );
		}
		else {
			switch(cov_.qp_solver_type_) {
				case QP_QPSOL:
				case QP_QPOPT:
			    case QP_QPSCHUR: {
					if( cov_.quasi_newton_==QN_LBFGS ) {
						matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
							= static_cast<IterQuantMatrixWithOpCreator*>(
								new IterQuantMatrixWithOpCreatorContinuous<MatrixSymPosDefLBFGS>()
								);
					}
					else if( cov_.quasi_newton_ == QN_PBFGS || cov_.quasi_newton_ == QN_LPBFGS ) {
						matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
							= static_cast<IterQuantMatrixWithOpCreator*>( 
								new IterQuantMatrixWithOpCreatorContinuous<
								ConstrainedOptimizationPack::MatrixHessianSuperBasicInitDiagonal>()
								);
					}
					else {
						matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
							= static_cast<IterQuantMatrixWithOpCreator*>( 
								new IterQuantMatrixWithOpCreatorContinuous<
								MatrixSymPosDefCholFactor>()
								);
					}
					break;
				}
				case QP_QPKWIK:
					matrix_iqa_creators[rSQPStateContinuousStorageMatrixWithOpCreator::Q_rHL]
						= new IterQuantMatrixWithOpCreatorContinuous<
									MatrixSymPosDefInvCholFactor>();
					break;
				default:
					assert(0);
			}
		}
				
		// Set the number of storage locations (these can be changed later)
		typedef rSQPStateContinuousStorage state_t;
		size_type storage_num[state_t::num_quantities];
		state_t::set_default_storage_num(storage_num);

		storage_num[state_t::Q_x]				= 2;
		storage_num[state_t::Q_f]				= 2;
		storage_num[state_t::Q_c]				= 2;

		storage_num[state_t::Q_opt_kkt_err]		= 2;
		storage_num[state_t::Q_feas_kkt_err]	= 2;
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

		// Now setup the reduced Hessian (this is kind of hacked but oh well)
		//
		// Here we will setup the matrix by setting it for the k-1 iteration
		// and then unsetting it.  This is a little bit of a hack but it should work.
		if( cov_.quasi_newton_ == QN_LBFGS ) {
			// Setup the LBFGS matrix.
			MatrixSymPosDefLBFGS *_rHL
				= dynamic_cast<MatrixSymPosDefLBFGS*>(&algo->rsqp_state().rHL().set_k(-1));
			if(_rHL) {
				_rHL->initial_setup(
					dof, cov_.num_lbfgs_updates_stored_
					,true  // Maintain the original matrix for QPOPT, QPSOL and QPSchur (iterative refinement)
					,cov_.qp_solver_type_==QP_QPSCHUR?true:false  // Maintian the inverse for QPSchur
					                                              // but now QPOPT or QPSOL
					,cov_.lbfgs_auto_scaling_
					);
				algo->rsqp_state().rHL().set_not_updated(-1);
			}
		}
		else if( cov_.quasi_newton_ == QN_PBFGS || cov_.quasi_newton_ == QN_LPBFGS ) {
			// Set the storage for the projected hessian for the
			// super basic variables B_RR.
			using ResourceManagementPack::ReleaseResource_ref_count_ptr;
			ConstrainedOptimizationPack::MatrixHessianSuperBasicInitDiagonal
				*_rHL = dynamic_cast<ConstrainedOptimizationPack::MatrixHessianSuperBasicInitDiagonal*>(
					&algo->rsqp_state().rHL().set_k(-1) );
			assert(_rHL); // Should not happen?
			ref_count_ptr<const MatrixSymWithOpFactorized>
				B_RR_ptr = NULL;
			if( cov_.quasi_newton_ == QN_LPBFGS ) {
				// Use a limited memory BFGS matrix initially
				ref_count_ptr<MatrixSymPosDefLBFGS>
					LB_RR_ptr = new MatrixSymPosDefLBFGS;
				LB_RR_ptr->initial_setup(
					dof, cov_.num_lbfgs_updates_stored_
					,true  // Maintain the original matrix for QPOPT, QPSOL and QPSchur (iterative refinement)
					,cov_.qp_solver_type_==QP_QPSCHUR?true:false  // Maintian the inverse for QPSchur
					                                              // but now QPOPT or QPSOL
					,cov_.lbfgs_auto_scaling_
					);
				LB_RR_ptr->init_identity( n-r, 1.0 ); // Must be sized when rHL is initialized
				B_RR_ptr = rcp::rcp_static_cast<const MatrixSymWithOpFactorized>(LB_RR_ptr);
			}
			else {
				// Go right ahead and use the dense projected BFGS matrix
				ref_count_ptr<MatrixSymPosDefCholFactor>
					DB_RR_ptr = new MatrixSymPosDefCholFactor(
						NULL    // Let it allocate its own memory
						,NULL,0 // ...
						,true   // Maintain the original matrix for QPOPT, QPSOL and QPSchur (iterative refinement)
						,cov_.qp_solver_type_==QP_QPSCHUR?true:false  // Maintian the cholesky factor for QPSchur
						                                              // but now QPOPT or QPSOL
						);
				DB_RR_ptr->init_identity( n-r, 1.0 ); // Must be sized when rHL is initialized
				B_RR_ptr = rcp::rcp_static_cast<const MatrixSymWithOpFactorized>(DB_RR_ptr);
			}
			// Finally initial the rHL matrix
			_rHL->initialize(
				dof, dof, NULL, NULL, NULL                     // n, n_R, i_x_free, i_x_fixed, bnd_fixed
				,B_RR_ptr                                      // B_RR_ptr (given ownership to destroy)
				,NULL, BLAS_Cpp::no_trans                      // B_RX_ptr, B_RX_trans
				,new SparseLinAlgPack::MatrixSymDiagonalStd()  // B_XX_ptr (sized to 0x0)
				);
			algo->rsqp_state().rHL().set_not_updated(-1);      // Set non updated so that it will be computed.
		}
		else if( cov_.quasi_newton_ == QN_BFGS ) {
			// Setup the dense, no frills, BFGS matrix for the cholesky factor
			switch( cov_.qp_solver_type_ ) {
			    case QP_QPSCHUR:
			    case QP_QPOPT:
			    case QP_QPSOL:
			    {
					MatrixSymPosDefCholFactor
						*_rHL = dynamic_cast<MatrixSymPosDefCholFactor*>(
							&algo->rsqp_state().rHL().set_k(-1) );
					assert(_rHL); // Should not happen?
					_rHL->init_setup(
						NULL    // Let it allocate its own memory
						,NULL,0 // ...
						,true // Maintain the original matrix for QPOPT, QPSOL and QPSchur (iterative refinement)
						,cov_.qp_solver_type_==QP_QPSCHUR?true:false // Maintian the cholesky factor for QPSchur
						                                             // but now QPOPT or QPSOL
					);
					
				}
			}
		}
		else {
			assert(0); // Local programming error?
		}
	}

	// /////////////////////////////////////////////////////
	// C.2. Create and set the decomposition system object

	DecompositionSystemVarReductStd*
		decomp_sys = NULL;

	if( !tailored_approach ) {
	
		// Only need a decomposition system object if we are not using the tailored approach

		if(trase_out)
			*trase_out << "\nCreating and setting the decomposition system object ...\n";

		if( !basis_sys_ptr_.get() ) {
			SparseCOOSolverCreator
				*sparse_solver_creator = NULL;
			if( cov_.direct_linear_solver_type_ == LA_MA28 ) {
				sparse_solver_creator
					= new SparseSolverPack::MA28SparseCOOSolverCreator(
							new SparseSolverPack::MA28SparseCOOSolverSetOptions
							, const_cast<OptionsFromStreamPack::OptionsFromStream*>(options_)
						);
			}
			else if ( cov_.direct_linear_solver_type_ == LA_MA48 ) {
				sparse_solver_creator
					= new SparseSolverPack::MA48SparseCOOSolverCreator(
							new SparseSolverPack::MA48SparseCOOSolverSetOptions
							, const_cast<OptionsFromStreamPack::OptionsFromStream*>(options_)
						);
			}
			else {
				assert(0);
			}

			// Note, the switch statement based code for the above if statements would not
			// compile under MipsPro 7.3.1.1m.

			basis_sys_ptr_ = new COOBasisSystem(sparse_solver_creator, true);
		}

		DecompositionSystemCoordinate* decomp_sys_aggr = 0;
		
		if ( cov_.null_space_matrix_type_ == NULL_SPACE_MATRIX_EXPLICIT ) {
			decomp_sys_aggr = new DecompositionSystemCoordinateDirect(
									basis_sys_ptr_.release(), true);
		}
		else if ( cov_.null_space_matrix_type_ == NULL_SPACE_MATRIX_IMPLICIT ) {
			decomp_sys_aggr = new DecompositionSystemCoordinateAdjoint(
									basis_sys_ptr_.release(), true);
		}
		else {
			assert(0);
		}

		assert(decomp_sys_aggr);
		typedef SparseSolverPack::TestingPack::BasisSystemTester basis_sys_tester_t;
		basis_sys_tester_t
			*basis_sys_tester = new basis_sys_tester_t;
		basis_sys_tester->solve_error_tol(1e-5); // ToDo: Make better documented?
		SparseSolverPack::TestingPack::BasisSystemTesterSetOptions
			basis_sys_options_setter( basis_sys_tester );
		basis_sys_options_setter.set_options( *options_ );
		decomp_sys_aggr->set_basis_sys_tester( basis_sys_tester ); // gives control to delete!

		// Note, the switch based code for the above if statements would not
		// compile under MipsPro 7.3.1.1m.

	   	decomp_sys = new DecompositionSystemVarReductStd( false, algo.get(), decomp_sys_aggr );

		algo->rsqp_state().set_decomp_sys(decomp_sys);

	}

	// /////////////////////////////////////////////////////
	// C.3  Create and set the step objects

	if(trase_out)
		*trase_out << "\nCreating and setting the step objects ...\n";

	{
		typedef ref_count_ptr<MeritFuncNLP> merit_func_ptr_t;
		merit_func_ptr_t merit_func;

		int step_num = 0;

		// Create the variable bounds testing object.  This will not
		// be used for algorithms where there are no variable bounds
		// but who cares.

		const value_type var_bounds_warning_tol = 1e-10;
		typedef ConstrainedOptimizationPack::VariableBoundsTester
			VariableBoundsTester;
		typedef rcp::ref_count_ptr<VariableBoundsTester> bounds_tester_ptr_t;
		bounds_tester_ptr_t
			bounds_tester = new VariableBoundsTester(
				  var_bounds_warning_tol			// default warning tolerance
				, algo_cntr.max_var_bounds_viol()	// default warning tolerance
				);
		{
			ConstrainedOptimizationPack::VariableBoundsTesterSetOptions
				options_setter( bounds_tester.get() );
			options_setter.set_options(*options_);
		}

		// Set the standard set of step objects.
		// This default algorithm is for NLPs with no variable bounds.

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
			if(options_) opt_setter.set_options( *options_ );

			algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
		}

		// (6) ReducedHessian
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

		// (6.-1) CheckSkipBFGSUpdate
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

		// (7) DepDirec
		if( !tailored_approach ) {
			algo->insert_step( ++step_num, DepDirec_name, new  DepDirecStd_Step );
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
					if(options_) direct_ls_opt_setter.set_options( *options_ );
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
					, new NewBasisSelectionStd_Strategy(decomp_sys)
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
	
		// Reconfigure the steps if the NLP has bounds

		if( algo_cntr.nlp().has_bounds() ) {

			if(trase_out)
				*trase_out << "\nReconfiguring steps for NLP with bounds ...\n";
			
			// If the NLP has bounds on the variables then we must replace
			// IndepDirec and move the convergence check to the end (along
			// with the calculation of the reduced gradient of the lagrangian).
			
			// Setup IndepDirec step
			if( cov_.qp_solver_type_ == QP_QPSCHUR || cov_.qp_solver_type_ == QP_QPKWIK ) {
			
				// Add iteration quantity for d_bounds
				algo->state().set_iter_quant( d_bounds_name
					, new IterQuantityAccessContinuous<SparseBounds>( 1, d_bounds_name ) );

				// Create the QP solver

				// RAB 8/28/00: In the future, all the QP solvers will also use this interface.
				typedef ref_count_ptr<QPSolverRelaxed> qp_solver_ptr_t;
				qp_solver_ptr_t qp_solver;

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
					default:
						assert(0);
				}

				// Create the QP solver tester object and set its options

				typedef ConstrainedOptimizationPack::QPSolverRelaxedTester qp_tester_t;
				typedef ref_count_ptr<qp_tester_t> qp_tester_ptr_t;

				qp_tester_ptr_t
					qp_tester = new qp_tester_t();
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

				Algorithm::poss_type poss;
				poss = algo->get_step_poss(IndepDirec_name);
				algo->remove_step( poss );	// remove any pre or post steps also
				algo->insert_step(
					  poss
					, IndepDirec_name
					, indep_direct_step	// Will clean up memory!
					);
				algo->insert_assoc_step(
					 poss
					, GeneralIterationPack::PRE_STEP
					, 1
					, "SetDBounds"
					, new  SetDBoundsStd_AddedStep
					);

			}
			else {

				// RAB 8/28/00: This step object and QP solver interface will be phased out!

				typedef ref_count_ptr<ReducedQPSolver> qp_solver_ptr_t;
				qp_solver_ptr_t qp_solver;

				switch(cov_.qp_solver_type_) {
					case QP_QPOPT:
					case QP_QPSOL:
					{
						ReducedQPSolverQPOPTSOL*  _qp_solver = 0;
						if(cov_.qp_solver_type_ == QP_QPOPT)
							_qp_solver = new ReducedQPSolverQPOPT;
						else
							_qp_solver = new ReducedQPSolverQPSOL;
						
						// If we are going to dump all of the iteration quantites we might as
						// well just print out the mappings to QPOPT and QPSOL.
						if(algo_cntr.journal_output_level() == PRINT_ITERATION_QUANTITIES)
							mapped_qp_file_ = mapped_qp_file_ptr_t( new std::ofstream("mapped_qp.txt") );
						_qp_solver->qp_mapping_output( mapped_qp_file_.get() );
						qp_solver = qp_solver_ptr_t( new ReducedQPSolverQPOPTSOLStd(_qp_solver,algo.get()) );
						break;
					}
				    default:
						assert(0);
				}

				Algorithm::poss_type poss;

				poss = algo->get_step_poss(IndepDirec_name);
				algo->replace_step(
							poss
							,new  IndepDirecExact_Step( 
								new ReducedQPSolverCheckOptimality(
									qp_solver, algo.get(), algo->algo_cntr().check_results() )
							  )
					    );
			}


			Algorithm::step_ptr_t calc_rgrad_lagr, calc_lambda, check_conv;

			// Remove and save CalcReducedGradLagr..., CalcLambdaIndep... and CheckConv...

			Algorithm::poss_type poss;
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

	}

	// 8/30/99: Add the iteration quantities for the QPSolverStats and 
	// ActSetStats and the added step that will calculate the
	// active set changes.

	{
		// Add active set and QP statistics to state object
		algo->state().set_iter_quant( act_set_stats_name
			, new IterQuantityAccessContinuous<ActSetStats>( 1, act_set_stats_name ) );
		algo->state().set_iter_quant( qp_solver_stats_name
			, new IterQuantityAccessContinuous<QPSolverStats>( 1, qp_solver_stats_name ) );

		if( algo_cntr.nlp().has_bounds() ) {

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
				new NewBasisSelectionStd_Strategy(decomp_sys) );
		CheckBasisFromPy_Step
			*basis_check_step2 = new CheckBasisFromPy_Step(
				new NewBasisSelectionStd_Strategy(decomp_sys) );
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

	// 7/3/00: Adding some steps to correct a bad initial guess.
	if( cov_.correct_bad_init_guess_ ) {
		typedef rcp::ref_count_ptr<CorrectBadInitGuessStd_AddedStep>
			corr_xinit_ptr_t;
		corr_xinit_ptr_t
			corr_xinit_ptr = new CorrectBadInitGuessStd_AddedStep;
		CorrectBadInitGuessStd_AddedStepSetOptions
			options_setter(corr_xinit_ptr.get());
		if(options_) options_setter.set_options(*options_);

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

}

void rSQPAlgo_ConfigMamaJama::init_algo(rSQPAlgoInterface& _algo)
{
	using DynamicCastHelperPack::dyn_cast;

	namespace rcp = ReferenceCountingPack;

	rSQPAlgo
#ifdef _WINDOWS
		&algo	= dynamic_cast<rSQPAlgo&>(_algo);
#else
		&algo	= dyn_cast<rSQPAlgo>(_algo);
#endif
	rSQPState	&state	= algo.rsqp_state();
	NLP			&nlp = algo.nlp();

	algo.max_iter( algo.algo_cntr().max_iter() );
	algo.max_run_time( algo.algo_cntr().max_run_time() );

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
	state.var_dep()			= Range1D(1, r);
	state.var_indep()		= Range1D(r + 1, n);
	state.con_decomp()		= Range1D(1, r);
	state.con_undecomp()	= (r < m) ? Range1D(r + 1, m) : Range1D(m + 1, m + 1);

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
	  const OptionsFromStreamPack::OptionsFromStream& options
	, SOptionValues *ov
	, std::ostream* trase_out
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
		const int num_opts = 14;
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
						ov->direct_linear_solver_type_ = LA_MA28;
					else if( linear_solver == "MA48" )
						ov->direct_linear_solver_type_ = LA_MA48;
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
						ov->qp_solver_type_ = QP_QPOPT;
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
	const SOptionValues& uov
	, SOptionValues *cov
	, std::ostream* trase_out )
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
				<< "\nqp_solver_type == AUTO: setting qp_solver_type = QPKWIK\n";
		cov->qp_solver_type_ = QP_QPKWIK;
	}
	else if(cov->qp_solver_type_ == QP_AUTO) {
		cov->qp_solver_type_ = uov.qp_solver_type_;
	}
	if( cov->line_search_method_ == LINE_SEARCH_AUTO && uov.line_search_method_ == LINE_SEARCH_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nline_search_method == AUTO: setting line_search_method = DIRECT\n";
		cov->line_search_method_ = LINE_SEARCH_DIRECT;
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
