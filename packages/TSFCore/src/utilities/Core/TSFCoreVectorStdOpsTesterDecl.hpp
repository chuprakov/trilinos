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

#ifndef TSFCORE_VECTOR_STD_OPS_TESTER_DECL_HPP
#define TSFCORE_VECTOR_STD_OPS_TESTER_DECL_HPP

#include "TSFCoreTypes.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace TSFCore {

/** \brief Testing class that tests all of the standard vector
 * operations defined in ??? using an arbitrary vector space.
 *
 * ToDo: Finish documenation!
 */
template <class Scalar>
class VectorStdOpsTester {

	///
	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

	/// Set the maximum relative error before an error is generated
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, maxRelErr );

	///
	VectorStdOpsTester(
		const ScalarMag   &maxRelErr = 100 * Teuchos::ScalarTraits<ScalarMag>::epsilon()
		);

	/** \brief Run the tests using a vector space.
	 *
	 * @param  vecSpc   [in] Vector space used to generate vectors in tests.
	 * @param  out      [in/out] If <tt>out!=NULL</tt> then <tt>*out</tt> will
	 *                  receive output about the tests.
	 * @param  dumpAll  [in] If <tt>true</tt> then vector elements will be printed after
	 *                  each transformation operation.  Default is <tt>false</tt>.
	 *
	 * @return Returns <tt>true</tt> if all of the tests check out and
	 * <tt>false</tt> otherwise.
	 */
	bool checkStdOps(
		const VectorSpace<Scalar>    &vecSpc
		,std::ostream                *out      = 0
		,const bool                  &dumpAll  = false
		);

};

} // namespace TSFCore

#endif // TSFCORE_VECTOR_STD_OPS_TESTER_DECL_HPP
