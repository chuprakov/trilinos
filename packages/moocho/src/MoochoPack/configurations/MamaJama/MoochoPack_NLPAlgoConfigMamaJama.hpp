// /////////////////////////////////////////////////////////////////////////
// rSQPAlgo_ConfigMamaJama.h
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
	  * group is "rSQPAlgo_ConfigMamaJama".  See the file rSQPpp.opt.rsqp_mama_jama_solve
	  * for these options and how they are specified.
	  *
	  * @param	options	[in] If NULL then no options will be set.  If !=NULL then this is the
	  * 				OptionsFromStream object that will be used to extract the
	  * 				options to use for the algorithm.  The state of this object must
	  * 				be maintained by the client until config_algo_cntr(...) is called
	  * 				and it is at this point that the options are read.
	  */
	void set_options( const OptionsFromStreamPack::OptionsFromStream* options );

	// ///////////////////
	// Overridden members

	///
	void config_algo_cntr(rSQPAlgoContainer& algo_cntr, std::ostream* trase_out);
	///
	void init_algo(rSQPAlgoInterface& algo);

public:

	// ///////////////////////////////////////////////////////////////////////
	// These should be private types but MS VC++ 6.0 will not allow it!

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
	enum EDirectLinearSolverType {
		LA_AUTO, LA_MA28, LA_MA48, LA_SUPERLU };
	///
	enum ENullSpaceMatrixType {
		NULL_SPACE_MATRIX_AUTO, NULL_SPACE_MATRIX_EXPLICIT
		, NULL_SPACE_MATRIX_IMPLICIT };
	///
	enum ERangeSpaceMatrixType {
		RANGE_SPACE_MATRIX_AUTO, RANGE_SPACE_MATRIX_COORDINATE
		, RANGE_SPACE_MATRIX_ORTHOGONAL };
	///
	enum EQuasiNewton {
		QN_AUTO, QN_BFGS, QN_PBFGS, QN_LBFGS, QN_LPBFGS };
	///
	enum EHessianInitialization {
		INIT_HESS_AUTO, INIT_HESS_IDENTITY, INIT_HESS_FIN_DIFF_SCALE_IDENTITY
		, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS };
	///
	enum EQPSolverType {
		QP_AUTO, QP_QPSOL, QP_QPOPT, QP_QPKWIK, QP_QPSCHUR };
	///
	enum ELineSearchMethod {
		LINE_SEARCH_AUTO, LINE_SEARCH_NONE, LINE_SEARCH_DIRECT
		, LINE_SEARCH_2ND_ORDER_CORRECT, LINE_SEARCH_WATCHDOG };
	///
	enum EMeritFunctionType {
		MERIT_FUNC_AUTO, MERIT_FUNC_L1, MERIT_FUNC_MOD_L1
		, MERIT_FUNC_MOD_L1_INCR };
	///
	enum EL1PenaltyParamUpdate {
		L1_PENALTY_PARAM_AUTO, L1_PENALTY_PARAM_WITH_MULT
		, L1_PENALTY_PARAM_MULT_FREE };

private:

	///
	struct SOptionValues {
		// Constructor (sets default values)
		SOptionValues();
		// Direct linear solvers
		EDirectLinearSolverType	direct_linear_solver_type_;
		// Variable Reduction,  Range/Null space decompositions
		ENullSpaceMatrixType	null_space_matrix_type_;
		ERangeSpaceMatrixType	range_space_matrix_type_;
		value_type				max_basis_cond_change_frac_;	// If < , don't change default
		// Reduced Hessian Approximations
		bool					exact_reduced_hessian_;
		EQuasiNewton			quasi_newton_;
		int						max_dof_quasi_newton_dense_;    // If < 0, don't change default
		int						num_lbfgs_updates_stored_;      // If < 0, don't change default
		bool					lbfgs_auto_scaling_;
		EHessianInitialization	hessian_initialization_;
		// QP subproblem solvers
		EQPSolverType			qp_solver_type_;
		bool                    reinit_hessian_on_qp_fail_;
		// Line search methods
		ELineSearchMethod		line_search_method_;
		EMeritFunctionType		merit_function_type_;
		EL1PenaltyParamUpdate	l1_penalty_param_update_;
		int						full_steps_after_k_;			// If < 0, do not use this option at all.
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
	SOptionValues		uov_;	// options set by user
	SOptionValues		cov_;	// current option values actually used

	// Keep a file that is used to output the mapped from a QP
	mapped_qp_file_ptr_t	mapped_qp_file_;

	// ///////////////////////////////////////////////////////
	// Private member functions

	/// Read in the options from a stream
	static void readin_options(
		const OptionsFromStreamPack::OptionsFromStream& options
		, SOptionValues *option_values, std::ostream* trase_out );

	/// Set the defaults for options not set by the user
	static void set_default_options(
		const SOptionValues& user_option_values
		, SOptionValues *current_option_values
		, std::ostream* trase_out );

};	// end class rSQPAlgo_ConfigMamaJama

}	// end namespace ReducedSpaceSQPPack 

#endif	// RSQP_ALGO_CONFIG_MAMA_JAMA_H
