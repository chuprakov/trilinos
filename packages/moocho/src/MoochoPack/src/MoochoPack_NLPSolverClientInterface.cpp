#include "../include/rSQPSolverClientInterface.h"

ReducedSpaceSQPPack::rSQPSolverClientInterface::rSQPSolverClientInterface(
		  int					max_iter				= 100
		, double				max_run_time			= 1e+10 // run forever
		, value_type			opt_tol					= 1e-6
		, value_type			feas_tol				= 1e-6
		, value_type			step_tol				= 1e-6
		, value_type			max_var_bounds_viol		= 0.0
		, EJournalOutputLevel 	journal_output_level 	= PRINT_NOTHING
		, int					journal_print_digits	= 6
		, bool					check_results			= true
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