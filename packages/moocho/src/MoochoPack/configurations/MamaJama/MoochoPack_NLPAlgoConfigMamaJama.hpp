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
  * It is assumed that there are explicit instantiations for the types it uses.
  */
class rSQPAlgo_ConfigMamaJama : public rSQPAlgo_Config {
public:

	/** @name Public Types */
	//@{

//	///
//  typedef	InvalidNLPType			InvalidNLPType;

	///
	enum EMeritFunctionType { MERIT_FUNC_L1, MERIT_FUNC_MOD_L1, MERIT_FUNC_MOD_L1_INCR };

	///
	enum ELineSearchMethod { LINE_SEARCH_NONE, LINE_SEARCH_DIRECT
		, LINE_SEARCH_2ND_ORDER_CORRECT, LINE_SEARCH_WATCHDOG };

	///
	enum EQPSolverType { QPSOL, QPOPT, QPKWIK, VE09, QPSCPD };

	///
	enum EFactorizationType { DIRECT_FACT, ADJOINT_FACT, AUTO_FACT };

	///
	enum EQuasiNewton { QN_AUTO, QN_BFGS, QN_LBFGS };

	///
	enum EHessianInitialization { INIT_HESS_IDENTITY, INIT_HESS_FIN_DIFF_SCALE_IDENTITY
		, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS };

	//@}


	/// Initialize

	///
	rSQPAlgo_ConfigMamaJama(
		  ReferenceCountingPack::ref_count_ptr<BasisSystem> basis_sys_ptr = 0
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			Gc_iq_creator_ptr = 0
		, ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
			U_iq_creator_ptr = 0
		);

	///
	~rSQPAlgo_ConfigMamaJama();

	/// Set the QPSolver type
	void qp_solver_type( EQPSolverType qp_solver_type )
	{	qp_solver_type_ = qp_solver_type; }

	/// Set the Factorzation type
	void factorization_type( EFactorizationType factorization_type )
	{	factorization_type_ = factorization_type; }

	/// Set whether extra tests are activated for checking intermediate results.
	void check_results( bool check_results )
	{	check_results_ = check_results; }

	/// Set whether we print the error in the QP subprolem..
	void print_qp_error( bool print_qp_error )
	{	print_qp_error_ = print_qp_error; }

	///
	/** Set the value of "big M" used in the solution of the relaxed QP for
	  * infeasible QPs.
	  */
	void bigM(value_type bigM)
	{	bigM_ = bigM;	}

	///
	/** Set the maximum change in the estimate of the condition number
	  * of the basis matrix before a change of basis is forced.
	  */
	void max_basis_cond_change_frac(value_type max_basis_cond_change_frac)
	{	max_basis_cond_change_frac_ = max_basis_cond_change_frac; }

	///
	/** Minumum fraction (0...1] of changes to active set before warm starts
	  * are used for QP subproblems with some QP solvers.
	  */
	void warm_start_frac(value_type warm_start_frac);

	///
	/** Set the type of merit function to use.
	  */
	void merit_function_type(EMeritFunctionType merit_function_type)
	{	merit_function_type_ = merit_function_type; }

	///
	/** Set the line search method to use.
	  */
	void line_search_method(ELineSearchMethod line_search_method)
	{	line_search_method_ = line_search_method; }

	///
	/** KKT error before the second order correction is attempted.
	  */
	void use_line_search_correct_kkt_tol(value_type use_line_search_correct_kkt_tol)
	{	use_line_search_correct_kkt_tol_ = use_line_search_correct_kkt_tol; }

	///
	/** Use full rSQP steps (alpha = 1) after full_steps_after_k rSQP iterations.
	  */
	void full_steps_after_k(int full_steps_after_k)
	{	full_steps_after_k_ = full_steps_after_k; }

	///
	/** Set the Quasi-Newton method to use for the reduced hessian approximation
	  *
	  * Here are the options:
	  * \begin{description}
	  *	\item[QN_AUTO] Limited memory is used if n-m > max_dof_quasi_newton_dense.
	  *	\item[QN_BFGS] Use the most appropriate dense form of the reduced hessian.
	  *	\iten[QN_LBFGS] Use the limited memory BFGS updating.
	  * \end{description}
	  *
	  * The default is QN_AUTO.
	  */
	void quasi_newton( EQuasiNewton quasi_newton )
	{	quasi_newton_ = quasi_newton; }

	///
	/** Set the initialization to be used for the reduced hessian approximation
	  *
	  * Here are the options:
	  * \begin{description}
	  *	\item[INIT_HESS_IDENTITY] Initialize the hessian to the identity matrix.
	  *	\item[INIT_HESS_FIN_DIFF_SCALE_IDENTITY] Peform a finite difference alone
	  *		the null space of the constriants then scale the identity matrix
	  *		by the largest finite difference term.
	  *	\iten[INIT_HESS_FIN_DIFF_SCALE_DIAGONAL] Peform a finite difference alone
	  *		the null space of the constriants then set the reduced hessain diagonal
	  *		to a modified finite difference.
	  * \end{description}
	  *
	  * The default is QN_AUTO.
	  */
	void hessian_initialization( EHessianInitialization hessian_initialization )
	{	hessian_initialization_ = hessian_initialization; }

	/// Set whether extra dampening is done on the quasi-Newton updates.
	void quasi_newton_dampening( bool quasi_newton_dampening )
	{	quasi_newton_dampening_ = quasi_newton_dampening; }

	///
	/** Set the maximum number of degrees of freedom for using
	  * a dense quasi-Newton matrix before a limited memory version
	  * is used.
	  *
	  * The default of max_dof_quasi_newton_dense = 1000 is used.
	  */
	void max_dof_quasi_newton_dense( int max_dof_quasi_newton_dense )
	{	max_dof_quasi_newton_dense_ = max_dof_quasi_newton_dense; }

	///
	/** If using limited memory BFGS, how many updates to store.
	  */
	void num_lbfgs_updates_stored( int num_lbfgs_updates_stored )
	{	num_lbfgs_updates_stored_ = num_lbfgs_updates_stored; }

	///
	/** If using limited memory BFGS, use auto rescaling of matrix or not.
	  */
	void lbfgs_auto_scaling( bool lbfgs_auto_scaling )
	{	lbfgs_auto_scaling_ = lbfgs_auto_scaling; }

	/// Set the options from stream object before configuring algorithm.
	void set_options( const OptionsFromStreamPack::OptionsFromStream* options )
	{	options_ = options; }

	// ///////////////////
	// Overridden members

	///
	void config_algo_cntr(rSQPAlgoContainer& algo_cntr, std::ostream* trase_out);
	///
	void init_algo(rSQPAlgoInterface& algo);

protected:
	typedef std::auto_ptr<std::ofstream> mapped_qp_file_ptr_t;
	typedef ReferenceCountingPack::ref_count_ptr<BasisSystem> basis_sys_ptr_t;
	typedef ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
		Gc_iq_creator_ptr_t;
	typedef ReferenceCountingPack::ref_count_ptr<IterQuantMatrixWithOpCreator>
		U_iq_creator_ptr_t;
	
	basis_sys_ptr_t		basis_sys_ptr_;	// Basis system object (if null will be set)
	Gc_iq_creator_ptr_t	Gc_iq_creator_ptr_;	// IQA creator for Gc
	U_iq_creator_ptr_t	U_iq_creator_ptr_;	// IQA creator for U which is N

	EQPSolverType		qp_solver_type_;
	EFactorizationType	factorization_type_;	// choosen by user
	EFactorizationType	fact_type_;				// actual type used.
	bool				check_results_;		// If true we will check all the results we can.
	bool				print_qp_error_;	// If true we print the QP error to journal.

	value_type			bigM_;
	value_type			max_basis_cond_change_frac_;	// default = 1.0, don't change default.
	value_type			warm_start_frac_;	// default = -1.0, don't change default
	EMeritFunctionType	merit_function_type_;
	ELineSearchMethod	line_search_method_;
	value_type 			use_line_search_correct_kkt_tol_;	// default = -1.0, don't change default
	int					full_steps_after_k_;	// default = -1, do use this option at all.
	EQuasiNewton		quasi_newton_;
	EHessianInitialization
						hessian_initialization_;
	bool				quasi_newton_dampening_;
	int					num_lbfgs_updates_stored_;
	bool				lbfgs_auto_scaling_;
	int					max_dof_quasi_newton_dense_;

	const OptionsFromStreamPack::OptionsFromStream
						*options_;
	// Keep a file that is used to output the mapped
	// from a QP
	mapped_qp_file_ptr_t	mapped_qp_file_;


};	// end class rSQPAlgo_ConfigMamaJama

}	// end namespace ReducedSpaceSQPPack 

#endif	// RSQP_ALGO_CONFIG_MAMA_JAMA_H
