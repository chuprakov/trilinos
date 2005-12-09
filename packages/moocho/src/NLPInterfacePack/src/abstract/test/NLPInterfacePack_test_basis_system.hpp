// //////////////////////////////////////////////////////////////
// NLPInterfacePack_test_basis_system.hpp
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

#ifndef TEST_BASIS_SYSTEM_H
#define TEST_BASIS_SYSTEM_H

#include "NLPInterfacePack_Types.hpp"

namespace OptionsFromStreamPack {
	class OptionsFromStream;
}

namespace NLPInterfacePack {

///
/** Test a \c BasisSystem object given matrices from a compatible \c NLPFirstOrder object.
 *
 * ToDo: Finish documentation!
 */
bool test_basis_system(
 	NLPFirstOrder                               *nlp
	,BasisSystem                                *basis_sys
	,OptionsFromStreamPack::OptionsFromStream   *options
	,std::ostream                               *out
	);

} // end NLPInterfacePack

#endif //TEST_BASIS_SYSTEM_H
