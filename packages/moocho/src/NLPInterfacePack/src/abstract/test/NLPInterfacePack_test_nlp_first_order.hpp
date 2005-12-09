// ////////////////////////////////////////////////////////////////////////
// NLPInterfacePack_test_nlp_first_order.hpp
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

#ifndef TEST_NLP_FIRST_ORDER_INFO_H
#define TEST_NLP_FIRST_ORDER_INFO_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"

namespace OptionsFromStreamPack {
	class OptionsFromStream;
}

namespace NLPInterfacePack {

///
/** Test an <tt>NLPFirstOrder</tt> object.
 *
 * @param  nlp     [in/out] %NLP object being tested.
 * @param  options [in] If <tt>options != NULL</tt> then the options to use are extracted
 *                 from <tt>*options</tt>.  If <tt>options == NULL</tt> then a default set
 *                 of options will be used that will be appropriate for even the largest %NLP
 *                 (see below).
 * @param  out     [in/out] If <tt>out != NULL</tt> then output will be set to <tt>*out</tt>.
 *                 The amount of output sent to <tt>*out</tt> depends on the options selected.
 *                 If <tt>out == NULL</tt> then no output is produced.
 *
 * This function uses the testing classes <tt>\ref AbstractLinAlgPack::VectorSpaceTester "VectorSpaceTester"</tt>
 * <tt>\ref NLPInterfacePack::NLPTester "NLPTester"</tt> and
 * <tt>\ref NLPInterfacePack::NLPFirstOrderInfoTester "NLPFirstOrderInfoTester"</tt> to perform many through tests
 * of an input <tt>\ref NLPInterfacePack::NLPFirstOrder "NLPFirstOrder"</tt> object.
 * The vector spaces exposed by <tt>\ref NLPInterfacePack::NLP "NLP"</tt> are thoroughly tested by the <tt>VectorSpaceTester</tt>
 * class.
 *
 * The options groups "VectorSpaceTester" (see <tt>\ref AbstractLinAlgPack::VectorSpaceTesterSetOptions "VectorSpaceTesterSetOptions"</tt>),
 * "%NLPTester" (see <tt>\ref NLPInterfacePack::NLPTesterSetOptions "NLPTesterSetOptions"</tt>), "%CalcFiniteDiffProd"
 * (see <tt>\ref NLPInterfacePack::CalcFiniteDiffProdSetOptions "CalcFiniteDiffProdSetOptions"</tt>) and "%NLPFirstOrderInfoTester"
 * (see <tt>\ref NLPInterfacePack::NLPFirstOrderInfoTesterSetOptions "NLPFirstOrderInfoTesterSetOptions"</tt>) are looked for in
 * in <tt>*options</tt> (if <tt>options != NULL</tt>) order to extract options to use for this testing function and the other testing
 * objects.
 */
bool test_nlp_first_order(
	NLPFirstOrder                               *nlp
	,OptionsFromStreamPack::OptionsFromStream   *options
	,std::ostream                               *out
	);

} // end namespace NLPInterfacePack

#endif // TEST_NLP_FIRST_ORDER_INFO_H
