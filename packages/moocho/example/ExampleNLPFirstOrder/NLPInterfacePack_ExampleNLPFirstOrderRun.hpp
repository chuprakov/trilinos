// ////////////////////////////////////////////////////////////////////
// ExampleNLPFirstOrderInfoRun.h
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

#ifndef EXAMPLE_NLP_FIRST_ORDER_INFO_RUN_H
#define EXAMPLE_NLP_FIRST_ORDER_INFO_RUN_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"

namespace NLPInterfacePack {

/** \defgroup ExampleNLPFirstOrderInfoRun_grp Helper function for ExampleNLPFirstOrderInfo */
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
 * @param  out         [in/out] If != NULL then *out gets the output (see the
 *                     options file "ExampleNLPFirstOrderInfoRun.opt").
 * @param  eout        [in/out] If != NULL then *eout gets minimal summary output.
 *
 * @returns true if the tests were successful, returns false otherwise.
 *
 * This function will read the file "ExampleNLPFirstOrderInfoRun.opt" in the
 * current directory to get the options to use.  The following is an example
 * of this file.
 * \verbinclude ExampleNLPFirstOrderInfoRun.opt
 * See <tt>\ref NLPInterfacePack::test_nlp_first_order_direct "test_nlp_first_order_direct()"</tt>
 * for descriptions of the options in this file.
 */
bool ExampleNLPFirstOrderInfoRun(
	const VectorSpace&   vec_space
	,value_type          xo
	,bool                has_bounds
	,bool                dep_bounded
	,std::ostream*       out
	,std::ostream*       eout
	);

//@}

} // end namespace NLPInterfacePack

#endif // EXAMPLE_NLP_FIRST_ORDER_INFO_RUN_H


