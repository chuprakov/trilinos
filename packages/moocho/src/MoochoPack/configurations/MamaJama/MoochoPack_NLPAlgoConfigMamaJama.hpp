// /////////////////////////////////////////////////////////////////////////
// rSQPAlgo_ConfigMamaJama.h

#ifndef RSQP_ALGO_CONFIG_MAMA_JAMA_H
#define RSQP_ALGO_CONFIG_MAMA_JAMA_H

#include <memory>
#include <fstream>

#include "../../include/rSQPAlgo_Config.h"
#include "../../include/rSQPAlgo.h"
#include "Misc/include/OptionsFromStream.h"

namespace ReducedSpaceSQPPack {

///
/** This is a do all configuration class for rSQPAlgo.
  *
  * 
  */
class rSQPAlgo_ConfigMamaJama : public rSQPAlgo_Config {
public:

	///
	rSQPAlgo_ConfigMamaJama(
		  ReferenceCountingPack::ref_count_ptr<BasisSystem>
			basis_sys_ptr = 0
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			Gc_iq_creator_ptr = 0
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			U_iq_creator_ptr = 0
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			HL_iq_creator_ptr = 0
		);

	///
	~rSQPAlgo_ConfigMamaJama();

	///
	/** Set the OptionsFromStream object that will be used for specifying the exact options.
	  *
	  * There are a lot of options associated with this configuration.  The main options
	  * group is "rSQPAlgo_ConfigMamaJama".  This options group and the default option
	  * values are shown below.
	  * 
	\begin{verbatim}

	options_group rSQPAlgo_ConfigMamaJama {

	*** Algorithmic/Implementation Options

	*** Direct fortran compatible linear solvers for the basis of the Jacobian
	    direct_linear_solver = MA28;

	*** Variable Reduction range/null space decomposition
	    null_space_matrix = AUTO;
	    range_space_matrix = AUTO;
	    max_basis_cond_change_frac = 1000; *** (+dbl) 

	*** Reduced Hessian Approximations
	    exact_reduced_hessian = false; *** Only if second order information is avalible.
	    quasi_newton = AUTO;
	    max_dof_quasi_newton_dense = 1000;
	    num_lbfgs_updates_stored = 7;
	    lbfgs_auto_scaling = true;
	    hessian_initialization = FINITE_DIFF_DIAGONAL_ABS;

	*** QP subproblem solver
	    qp_solver = QPKWIK;
	    warm_start_frac = 1.0;

	*** Line search methods
	    line_search_method = DIRECT;
	    merit_function_type = L1;
	    L1_penalty_parameter_update = WITH_MULT;

	}
	\end{verbatim}
	  *
	  * These options are described here:
	  * 
	  * \begin{enumeration}
	  * \item {\bf Direct Linear Solvers : }  These are options that the determine the software
	  * 	used to perform a sparse direct factorization of the basis of the Jacobian (Gc') of
	  * 	the constraints.
	  * 	\begin{description}
	  * 	\item[MA28] The Harwell package MA28
	  * 	\item[MA48] The Harwell package MA48.  Note that this also requires MA28 in order to
	  * 		find a nonsingular basis since this is not directly supported in MA48.
	  * 	\end{description}
	  *	\item {\bf Variable Reduction Range/Null Space Decomposition : } These options determine
	  *		the specifics of how the range (Y) and null (Z) space matrices are represented.
	  *		\begin{description}
	  *		\item[null_space_matrix] Determines how the variable reduction matrix
	  *			#Z = [ D; I ]# is implemented where #Gc' = [ C , N ], D = -inv(C)*N#.
	  *			\begin{description}
	  *			\item[AUTO] Let the algorithm decide what to do (this is usually prefered).
	  *			\item[EXPLICIT] Store #D = -inv(C)*N# explicitly and then perform all computation with it.
	  *			\item[EXPLICIT] Implement all operations with #C# and #N# (less storage).
	  *			\end{description}
	  *		\item[ToDo: Finish This!]
	  *		\end{description}
	  *	\item {\bf ToDo: Finish This!}
	  * \end{enumeration}
	  *
	  * In addition to the above options group the following options groups can be included and they are:
	  *
	  * If a "Tailored Approach" NLP is being solved,
	  * the options group "EvalNewPointTailoredApproachStd" can be inclded
	  * (see the class \Ref{EvalNewPointTailoredApproachStd_StepSetOptions})
	  * to change the options for the evaluation of the NLP point for a "Tailored
	  * Approach" NLP.  Also, options for the finite difference tester used in this step
	  * can be set by including the options group "NLPrSQPTailoredApproachTester"
	  * (see \Ref{NLPrSQPTailoredApproachTesterSetOptions}).
	  *
	  * If quasi-Newton BFGS updating is being performed,
	  * the options group "InitFinDiffReducedHessian" can be included
	  * (see the class \Ref{InitFinDiffReducedHessian_StepSetOptions})
	  * to change options for the reduced hessian initialization.
	  * 
	  * If quasi-Newton BFGS updating is being performed,
	  * the options group "CheckSkipBFGSUpdateStd" can be included
	  * (see the class \Ref{CheckSkipBFGSUpdateStd_StepSetOptions})
	  * to change options for the skipping of the BFGS update.
	  *
	  * The options group "CheckConvergenceStd_AddedStep" can be
	  * included (see the class \Ref{CheckConvergenceStd_AddedStepSetOptions})
	  * to adjust some of the finer details of the convergence criteria.
	  *
	  * To change options for the direct line search for the SQP step
	  * (if a line seach method is being used)
	  * the options group "DirectLineSearchArmQuadSQPStep" can be included
	  * (see the class \Ref{DirectLineSearchArmQuad_StrategySetOptions}).
	  * This direct line search object is used by all of the line search.
	  * 
	  * If the watchdog line search method has been selected,
	  * the options group "LineSearchWatchDog" can be included
	  * (see the class \RefLineSearchWatchDog_AddedStepSetOptions})
	  * to change the options for the watchdog line search.
	  *
	  * If an L1 merit function is being used (several different types),
	  * the options group "MeritFuncPenaltyParamUpdate" can be included
	  * (see the class \Ref{MeritFunc_PenaltyParamUpdate_AddedStepSetOptions})
	  * to change the options for the updates of the penalty parameters.
	  *
	  * If the option "merit_function_type = MODIFIED_L1_INCR" is selected,
	  * the options group "MeritFuncModifiedL1LargerSteps" can be included
	  * (see the class \Ref{MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions})
	  * to change the options for the updates of the penalty parameters
	  * for the modified L1 merit function only.
	  *
	  * If the second order correction line search is being used
	  * (line_search_method = 2ND_ORDER_CORRECT) for the SQP step
	  * the options groups "LineSearch2ndOrderCorrect"
	  * (see the class \Ref{LineSearch2ndOrderCorrect_StepSetOptions}) and
	  * "DirectLineSearchArmQuad2ndOrderCorrectNewton"
	  * (see the class \Ref{DirectLineSearchArmQuad_StrategySetOptions})
	  * can be included.  See the class \Ref{LineSearch2ndOrderCorrect_Step}
	  * for a description of these.
	  *
	  * If the QP solver QPSCPD has been selected (qp_solver = QPSCPD)
	  * , the options group "QPSCPD" can be included
	  * (see the class \Ref{QPMixedFullReducedQPSCPDSolverSetOptions})
	  * to change some of the default options with this QP solver.
	  * 
	  * If the MA28 linear solver is selected (direct_linear_solver = MA28) then
	  * the options group "MA28SparseCOOSolver" can be included.
	  *
	  * @param	options	[I]	If NULL then no options will be set.  If !=NULL then this is the
	  * 					Options from stream object that will be used to extract the
	  * 					options to use for the algorithm.  The state of this object must
	  * 					be maintained by the client until config_algo_cntr(...) is called
	  * 					and it is at this point that the options are read.
	  */
	void set_options( const OptionsFromStreamPack::OptionsFromStream* options );

	// ///////////////////
	// Overridden members

	///
	void config_algo_cntr(rSQPAlgoContainer& algo_cntr, std::ostream* trase_out);
	///
	void init_algo(rSQPAlgoInterface& algo);

private:

	// ///////////////////////////////////////////////////////////////////////
	// Private types

	///
	typedef std::auto_ptr<std::ofstream> mapped_qp_file_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<BasisSystem> basis_sys_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
		Gc_iq_creator_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
		U_iq_creator_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			HL_iq_creator_ptr_t;
	///
	enum EMeritFunctionType { MERIT_FUNC_L1, MERIT_FUNC_MOD_L1, MERIT_FUNC_MOD_L1_INCR };
	///
	enum EL1PenaltyParamUpdate { L1_PENALTY_PARAM_WITH_MULT, L1_PENALTY_PARAM_MULT_FREE };
	///
	enum ELineSearchMethod { LINE_SEARCH_NONE, LINE_SEARCH_DIRECT
		, LINE_SEARCH_2ND_ORDER_CORRECT, LINE_SEARCH_WATCHDOG };
	///
	enum EQPSolverType { QPSOL, QPOPT, QPKWIK, QPSCPD };
	///
	enum EDirectLinearSolverType { MA28, MA48 };
	///
	enum ENullSpaceMatrixType { NULL_SPACE_MATRIX_EXPLICIT, NULL_SPACE_MATRIX_IMPLICIT
		, NULL_SPACE_MATRIX_AUTO };
	///
	enum ERangeSpaceMatrixType { RANGE_SPACE_MATRIX_COORDINATE, RANGE_SPACE_MATRIX_ORTHOGONAL };
	///
	enum EQuasiNewton { QN_AUTO, QN_BFGS, QN_LBFGS };
	///
	enum EHessianInitialization { INIT_HESS_IDENTITY, INIT_HESS_FIN_DIFF_SCALE_IDENTITY
		, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS };
	///
	struct SOptionValues {
		// constructor (default values)
		SOptionValues();
		// Direct linear solvers
		EDirectLinearSolverType	direct_linear_solver_type_;
		// Variable Reduction,  Range/Null space decompositions
		ENullSpaceMatrixType	null_space_matrix_type_;		// set by user
		ENullSpaceMatrixType	null_space_matrix_type_used_;	// actually used
		ERangeSpaceMatrixType	range_space_matrix_type_;
		value_type				max_basis_cond_change_frac_;	// default = -1.0, don't change default.
		// Reduced Hessian Approximations
		bool					exact_reduced_hessian_;			// default = false
		EQuasiNewton			quasi_newton_;					// set by user
		EQuasiNewton			quasi_newton_used_;				// actually used
		int						max_dof_quasi_newton_dense_;
		int						num_lbfgs_updates_stored_;
		bool					lbfgs_auto_scaling_;
		EHessianInitialization	hessian_initialization_;
		// QP subproblem solvers
		EQPSolverType			qp_solver_type_;
		value_type				warm_start_frac_;				// default = -1.0, don't change default
		// Line search methods
		ELineSearchMethod		line_search_method_;
		EMeritFunctionType		merit_function_type_;
		EL1PenaltyParamUpdate	l1_penalty_param_update_;
		int						full_steps_after_k_;			// default = -1, do use this option at all.
		// Not used for now
		bool					print_qp_error_;				// If true we print the QP error to journal.
		bool					correct_bad_init_guess_;
	};

	/// Possible user supplied stuff
	basis_sys_ptr_t		basis_sys_ptr_;	// Basis system object (if null will be set)
	Gc_iq_creator_ptr_t	Gc_iq_creator_ptr_;	// IQA creator for Gc
	U_iq_creator_ptr_t	U_iq_creator_ptr_;	// IQA creator for U which is N
	HL_iq_creator_ptr_t	HL_iq_creator_ptr_;	// IQA creator for HL

	/// Pointer to options
	const OptionsFromStreamPack::OptionsFromStream
						*options_;
	/// Options
	SOptionValues		cov_;	// current option values

	// Keep a file that is used to output the mapped
	// from a QP
	mapped_qp_file_ptr_t	mapped_qp_file_;

	// ///////////////////////////////////////////////////////
	// Private member functions

	/// Read in the options from a stream
	static void readin_options( const OptionsFromStreamPack::OptionsFromStream& options
		, SOptionValues *option_values, std::ostream* trase_out );

};	// end class rSQPAlgo_ConfigMamaJama

}	// end namespace ReducedSpaceSQPPack 

#endif	// RSQP_ALGO_CONFIG_MAMA_JAMA_H
