// ////////////////////////////////////////////////////////////////////
// ExampleNLPFirstOrderDirectRun.h
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
//

#ifndef EXAMPLE_NLP_FIRST_ORDER_DIRECT_RUN_H
#define EXAMPLE_NLP_FIRST_ORDER_DIRECT_RUN_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"
#include "ReducedSpaceSQPPack/Configurations/rSQPppSolver.h"

namespace NLPInterfacePack {

/** \defgroup ExampleNLPFirstOrderDirectRun_grp Helper function for ExampleNLPFirstOrderDirect */
//@{

///
/** Function accepts a VectorSpace object and then uses it to define
 * an example NLP and run a whole battery of tests.
 *
 * @param  vec_space   [in] The vector space object used to create all of the
 *                     needed vector spaces and vectors.  This vector space and
 *                     the vectors it creates will get a though testing.
 * @param  xo          [in] The initial starting point for unknown variables (before
 *                     they are forced in bounds).
 * @param  has_bounds  [in] If true, then the NLP will have bounds on the variables.
 * @param  dep_bounded [in] (valid only if has_bounds == true) If true, then
 *                     the dependent variables will be bounded, if false the
 *                     independent variables will be bounded.
 * @param  console_out [in/out] If != NULL then *console_out gets the output.
 * @param  error_out   [in/out] If != NULL then *eout gets minimal summary output.
 * @param  throw_solve_exception
 *                     [in] If true then solver will not throw exception (but other code may).
 * @param  algo_out    [in/out] If != NULL then it gets algo outptut, otherwise goes to 'rSQPppAlgo.out'
 * @param  summary_out [in/out] If != NULL then it gets summary outptut, otherwise goes to 'rSQPppSummary.out'
 * @param  journal_out [in/out] If != NULL then it gets journal outptut, otherwise goes to 'rSQPppJournal.out'
 *
 * @returns Returns the return value from <tt>ReducedSpaceSQPPack::rsqp_mama_jama_solve()</tt>
 * (see this function for most of the documentation).
 */
ReducedSpaceSQPPack::rSQPppSolver::ESolutionStatus
ExampleNLPFirstOrderDirectRun(
	const VectorSpace&   vec_space
	,value_type          xo
	,bool                has_bounds
	,bool                dep_bounded
	,std::ostream*       console_out
	,std::ostream*       error_out
	,bool                throw_sole_exception = false
	,std::ostream*       algo_out             = NULL
	,std::ostream*       summary_out          = NULL
	,std::ostream*       journal_out          = NULL
	);

//@}

} // end namespace NLPInterfacePack

#endif // EXAMPLE_NLP_FIRST_ORDER_DIRECT_RUN_H


