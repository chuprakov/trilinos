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

// ///////////////////////////////////////////////////////////////////
// TSFCoreLinearOpTesterDecl.hpp

#ifndef TSFCORE_LINEAR_OP_TESTER_DECL_HPP
#define TSFCORE_LINEAR_OP_TESTER_DECL_HPP

#include "TSFCoreTypes.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace TSFCore {

///
/** Testing class for <tt>LinearOp</tt>.
 *
 * This testing class performs a many tests as possible just given a
 * <tt>LinearOp</tt> object using the function <tt>check()</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_ANA_Development_grp
 */
template<class Scalar>
class LinearOpTester {
public:

	///
	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  ///
	/** Set the tolerance above which a relative error will generate a warning message.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, warning_tol )
  ///
	/** Set the tolerance above which a relative error will generate a error message and result in test failure.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, error_tol )

	///
	/** Default constructor.
	 */
	LinearOpTester(
		const ScalarMag      warning_tol = 1e-13
		,const ScalarMag     error_tol   = 1e-10
		);

	///
	/** Check out a linear operator.
	 *
	 * @param  op    [in] The linear operator to check out
	 * @param  out   [in/out] If <tt>out!=NULL</tt> then trace output
	 *               about the tests performed will be sent to <tt>*out</tt>.
	 *
	 * This function performs a number of tests on <tt>op</tt>:<ul>
	 *
	 * <li>Checks that the domain and range spaces are valid
	 *
	 * <li>Creates temporary vectors using the domain and range spaces
	 *
	 * <li>Checks that the non-transposed and the transposed operator agree.
	 *     The operator and adjoint operator must obey the defined
	 *     scalar product.  Specifically, for any two vectors \f$w\f$ (in the
	 *     domain space \f$\mathcal{D}\f$) and \f$u\f$ (in the range space
	 *     \f$\mathcal{R}\f$) the adjoint operation must obey:
	 *     \f[<u,A v>_{\mathcal{R}} = <A^T u, v>_{\mathcal{D}}\f]
   *     where \f$<.,.>_{\mathcal{R}}\f$ is the scalar product defined by
	 *     <tt>op.range()->scalarProd()</tt> and \f$<.,.>_{\mathcal{D}}\f$
	 *     is the scalar product defined by
	 *     <tt>op.domain()->scalarProd()</tt>.
	 *
	 * </ul>
	 *
	 * All relative errors that exceed <tt>warning_tol()</tt> but do not
	 * exceed <tt>error_tol</tt> will result in special warning messages
	 * printed to <tt>*out</tt> (if <tt>out!=NULL</tt>).
	 *
	 * @return The function returns <tt>true</tt> if all of the tests
	 * where within the <tt>error_tol()</tt> and returns <tt>false</tt>
	 * if not.
	 *
	 * The best way to see what this testing function is doing is to run
	 * the test with <tt>out!=NULL</tt> and to look at the
	 * implementation by clicking on the following link to the source code:
	 */
	bool check(
		const LinearOp<Scalar>      &op
		,std::ostream               *out
		) const;

}; // class LinearOpTester

} // namespace TSFCore

#endif // TSFCORE_LINEAR_OP_TESTER_DECL_HPP
