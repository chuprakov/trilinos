// /////////////////////////////////////////////////////////////////////////////
// rSQPSolverClientInterface.cpp
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

#include "ReducedSpaceSQPPack/include/rSQPSolverClientInterface.h"
#include "GeneralIterationPack/include/AlgorithmTrack.h"

ReducedSpaceSQPPack::rSQPSolverClientInterface::rSQPSolverClientInterface(
	int                      max_iter
	,double                  max_run_time
	,value_type              opt_tol
	,value_type              feas_tol
	,value_type              comp_tol
	,value_type              step_tol
	,EJournalOutputLevel     journal_output_level
	,int                     journal_print_digits
	,bool                    check_results
	)
	:max_iter_(max_iter)
	,max_run_time_(max_run_time)
	,opt_tol_(opt_tol)
	,feas_tol_(feas_tol)
	,comp_tol_(comp_tol)
	,step_tol_(step_tol)
	,journal_output_level_(journal_output_level)
	,journal_print_digits_(journal_print_digits)
	,check_results_(check_results)
{}
