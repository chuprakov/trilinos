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

#include "ReducedSpaceSQPPack/src/rSQPAlgo_Config.h"
#include "ReducedSpaceSQPPack/src/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/Configurations/shared/DecompositionSystemStateStepBuilderStd.h"
#include "OptionsFromStream.h"

namespace ReducedSpaceSQPPack {

///
/** This is a do all configuration class for <tt>rSQPAlgo</tt>.
 *
 * This class relies on the builder class <tt>DecompositionSystemStateStepBuilderStd</tt>
 * to perform many different tasks.
 *
 * Options specific for to this configuration class and the classes that
 * it works with that can be set through <tt>this->set_options()</tt>, see the files
 * <tt>\ref DecompositionSystemStateStepBuilderStd_opts "rSQPpp.opt.DecompositionSystemStateStepBuilderStd"</tt>.
 * and <tt>\ref rSQPAlgo_ConfigMamaJama_opts "rSQPpp.opt.rSQPAlgo_ConfigMamaJama"</tt>.
 *
 * Note that all built-in support for basis permutations and direct sparse solvers
 * can be left out if the macro RSQPPP_NO_BASIS_PERM_DIRECT_SOLVERS is defined
 * This will result in smaller executables for programs that don't need this
 * extra functionality..
 * 
 * ToDo: Finish documentation!
 */
class rSQPAlgo_ConfigMamaJama : public rSQPAlgo_Config {
public:

	///
	rSQPAlgo_ConfigMamaJama();

	///
	~rSQPAlgo_ConfigMamaJama();

	/** Overridden from rSQPAlgo_Config */
	//@{

	///
	/** Set the <tt>OptionsFromStream</tt> object that will be used for specifying the options.
	 *
	 *  @param  options
	 *               [in] If \c NULL then no options will be set.  If <tt>!=NULL</tt> then
	 *               this is the \c OptionsFromStream object that will be used to extract the
	 *               options to use for the algorithm.  The state of this object must
	 *               be maintained by the client until \c config_algo_cntr() is called
	 *               and it is at this point that the options are read.
	 *
	 */
	void set_options( const options_ptr_t& options );
	///
	const options_ptr_t& get_options() const;
	///
	void config_algo_cntr(rSQPAlgoContainer* algo_cntr, std::ostream* trase_out);
	///
	void init_algo(rSQPAlgoInterface* algo);

	//@}

public:

	/** @name Enums for variaous options categories */
	//@{

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
		, LINE_SEARCH_2ND_ORDER_CORRECT, LINE_SEARCH_WATCHDOG
		, LINE_SEARCH_FILTER };
	///
	enum EMeritFunctionType {
		MERIT_FUNC_AUTO, MERIT_FUNC_L1, MERIT_FUNC_MOD_L1
		, MERIT_FUNC_MOD_L1_INCR };
	///
	enum EL1PenaltyParamUpdate {
		L1_PENALTY_PARAM_AUTO, L1_PENALTY_PARAM_WITH_MULT
		, L1_PENALTY_PARAM_MULT_FREE };

	//@}

	/** @name Struct for options values */
	//@{

	///
	struct SOptionValues {
		// Constructor (sets default values)
		SOptionValues();
		// Variable Reduction,  Range/Null space decompositions
		value_type				max_basis_cond_change_frac_;	// If < , don't change default
		// Reduced Hessian Approximations
		bool					exact_reduced_hessian_;
		EQuasiNewton			quasi_newton_;
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
	};

	//@}

private:

	/// Builder class for some common code
	DecompositionSystemStateStepBuilderStd   decomp_sys_step_builder_;

	/// Smart pointer to options
	options_ptr_t      options_;

	/// Options structs
	SOptionValues       uov_; // options set by user
	SOptionValues       cov_; // current option values actually used

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

/** \defgroup rSQPAlgo_ConfigMamaJama_opts Options for rSQPAlgo_ConfigMamaJama.
 *
 * The following is the contents of the file <tt>rSQPpp.opt.rSQPAlgo_ConfigMamaJama</tt>
 * which are options specific to the class <tt>ReducedSpaceSQPPack::rSQPAlgo_ConfigMamaJama</tt>
 * and the class objects that it configures.
 *
 * \verbinclude rSQPpp.opt.rSQPAlgo_ConfigMamaJama
 */

}	// end namespace ReducedSpaceSQPPack 

#endif	// RSQP_ALGO_CONFIG_MAMA_JAMA_H
