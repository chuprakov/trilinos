// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// /////////////////////////////////////////////////////////////////////////
// TSFCoreTestingToolsDecl.hpp

#ifndef TSFCORE_TESTING_TOOLS_DECL_HPP
#define TSFCORE_TESTING_TOOLS_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

/** \defgroup TSFCore_test_tools_code_grp Miscellaneous C++ utility code for testing and debugging.

\ingroup TSFCore_ANA_Development_grp

Here is some assorted C++ code to aid in testing and debugging
%TSFCore code.

*/

///
/** Return "passed" or "failed".
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
inline
const char* passfail(bool pass) { return pass ? "passed" : "failed"; }

///
/** Return relative error of two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
template <class Scalar>
Scalar relErr( const Scalar &s1, const Scalar &s2 );

///
/** Compute, check and optionally print the relative error in two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
template<class Scalar>
bool testRelErr(
	const std::string    &v1_name
	,const Scalar        &v1
	,const std::string   &v2_name
	,const Scalar        &v2
	,const std::string   &maxRelErr_name
	,const Scalar        &maxRelErr
	,std::ostream        *out
	);

///
/** Output operator to pretty print any <tt>TSFCore::Vector</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const Vector<Scalar>& v );

///
/** Output operator to pretty print any <tt>TSFCore::LinearOp</tt> object.
 *
 * Calls <tt>operator<<( std::ostream, const LinOpNonPersisting<Scalar> )</tt>
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const LinearOp<Scalar>& M );

///
/** Output operator to pretty print any <tt>TSFCore::LinOpNonPersisting</tt> object.
 *
 * Calls <tt>operator<<( std::ostream, const LinOpNonPersisting<Scalar> )</tt>
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const LinOpPersisting<Scalar>& M );

///
/** Output operator to pretty print any <tt>TSFCore::LinOpNonPersisting</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_test_tools_code_grp
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const LinOpNonPersisting<Scalar>& M );

} // namespace TSFCore

#endif // TSFCORE_TESTING_TOOLS_DECL_HPP
