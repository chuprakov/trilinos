#include "../include/rSQPSolverClientInterface.h"

ReducedSpaceSQPPack::rSQPSolverClientInterface::rSQPSolverClientInterface(
		  int					max_iter
		, double				max_run_time
		, value_type			opt_tol
		, value_type			feas_tol
		, value_type			step_tol
		, value_type			max_var_bounds_viol
		, EJournalOutputLevel 	journal_output_level
		, int					journal_print_digits
		, bool					check_results
		)
	:
		  max_iter_(max_iter)
		, max_run_time_(max_run_time)
		, opt_tol_(opt_tol)
		, feas_tol_(feas_tol)
		, step_tol_(step_tol)
		, max_var_bounds_viol_(max_var_bounds_viol)
		, journal_output_level_(journal_output_level)
		, journal_print_digits_(journal_print_digits)
		, check_results_(check_results)
{}