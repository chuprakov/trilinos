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

// ///////////////////////////////////////////////////////////////
// TSFCoreSolversNormedConvergenceTesterDecl.hpp

#ifndef TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
#define TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP

#include "TSFCoreSolversAttachConvergenceTesterBase.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Convergence test based on a relative error tolerance.
 *
 * This concrete subclass is derived from
 * <tt>AttachedConvergenceTesterBase</tt> and therefore inherits all
 * of the machinary for handling attachec convergence tests.
 *
 * 
 */
template<class Scalar>
class NormedConvergenceTester : public AttachConvergenceTesterBase<Scalar> {
public:

	///
	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  ScalarMagnitude;
	///
	typedef typename AttachConvergenceTesterBase<Scalar>::EAttachmentMode EAttachmentMode;

	/** @name Constructors / initializers */
	//@{
	
	///
	/** Calls <tt>tol()</tt> and <tt>attachmentMode()</tt>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->tol() == tol</tt>
	 * <li><tt>this->attachmentMode()==attachmentMode</tt>.
	 * </ul>
	 *
	 * By default <tt>*this</tt> is constructed to ignore an attached
	 * convergence tester.
	 */
	NormedConvergenceTester(
		const ScalarMagnitude tol             = ScalarMagnitude(1e-12)
		,const EAttachmentMode attachmentMode = AttachConvergenceTesterBase<Scalar>::ATTACHED_TEST_EXCLUDE
		);

	///
	/** Set the tolerance used for the convergence check.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->tol() == tol</tt>
	 * <li><tt>this->minMaxErr() == ScalarMagnitude(1e+50)</tt>
	 * </ul>
	 */
	void tol(
		const ScalarMagnitude tol
		);
	
	/// Return the tolerance passed into <tt>initialize()</tt>.
	ScalarMagnitude tol() const;

	/// Return the minimum (over all iterations) maximum (over all right-hand sides) error seen.
	ScalarMagnitude minMaxError() const;

	/// Set the value returned by <tt>minMaxError()</tt> as a hack (should not be called generally)
	void minMaxErr( const ScalarMagnitude& minMaxErr );

	//@}

protected:

	/** @name Overridden from AttachedConvergenceTesterBase */
	//@{

	///
	void protectedReset();
	///
	void protectedConvStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);

	//@}

private:

	ScalarMagnitude                     tol_;
	ScalarMagnitude                     minMaxErr_;

	std::valarray<ScalarMagnitude>      norms_;          // cache

}; // class NormedConvergenceTester

// ///////////////////////////
// Inline members

template<class Scalar>
inline
typename NormedConvergenceTester<Scalar>::ScalarMagnitude
NormedConvergenceTester<Scalar>::tol() const
{
	return tol_;
}

template<class Scalar>
inline
typename NormedConvergenceTester<Scalar>::ScalarMagnitude
NormedConvergenceTester<Scalar>::minMaxError() const
{
	return minMaxErr_;
}

template<class Scalar>
inline
void NormedConvergenceTester<Scalar>::minMaxErr( const ScalarMagnitude& minMaxErr )
{
	minMaxErr_ = minMaxErr;
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
