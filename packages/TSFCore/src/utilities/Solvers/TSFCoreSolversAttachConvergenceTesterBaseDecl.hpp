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
// TSFCoreSolversAttachConvergenceTesterBaseDecl.hpp

#ifndef TSFCORE_SOLVERS_ATTACH_CONVERGENCE_TESTER_BASE_DECL_HPP
#define TSFCORE_SOLVERS_ATTACH_CONVERGENCE_TESTER_BASE_DECL_HPP

#include "TSFCoreSolversConvergenceTester.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Node base subclass for all <tt>ConvergenceTester</tt> subclasses
 * that just want default functionality for the attachment of
 * convergence testers.
 *
 * This subclass sets up the machinary for setting up an attached
 * convergence tester and considering it in testing.
 *
 * <b>Note to subclass developers:</b>
 *
 * The only functions that must be overridden are the public function
 * protected functions <tt>protectedReset()</tt> and
 * <tt>protectedConvStatus()</tt>.
 */
template<class Scalar>
class AttachConvergenceTesterBase : public ConvergenceTester<Scalar> {
public:

	/** @name Public types */
	//@{

	///
	enum EAttachmentMode {
		ATTACHED_TEST_EXCLUDE  ///< The attached convergence tester will never even be called.
		,ATTACHED_TEST_IGNORE  ///< The attached convergence tester will be called but the result ignored.
		,ATTACHED_TEST_INSTEAD ///< The attached convergence tester will be used in place of this convergence test.
		,ATTACHED_TEST_AND     ///< The attached convergence tester will be called and have its results and'ed.
		,ATTACHED_TEST_OR      ///< The attached convergence tester will be called and have its results or'ed
	};

	//@}

	///
	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  ScalarMagnitude;

	/** @name Constructors / initializers */
	//@{

	/// Calls <tt>attachmentMode()</tt>
	AttachConvergenceTesterBase(
		const EAttachmentMode attachmentMode = ATTACHED_TEST_AND
		);

	/// Set the attachment mdoe (see <tt>EAttachmentMode</tt>)
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EAttachmentMode, attachmentMode )

	//@}

	/** @name Overridden from ConvergenceTester */
	//@{

	///
	/** Overridden.
	 *
	 * Overridden to call <tt>protectedReset()</tt> and will call
	 * <tt>getAttachedConvTester()->reset()<tt> if
	 * <tt>getAttachedConvTester().get()!=NULL &&
	 * attachmentMode()!=ATTACH_TEST_EXCLUDE</tt>.
	 */
	void reset();
	///
	/** Overridden.
	 *
	 * Overridden to call <tt>protectedConvStatus()</tt> and will call
	 * <tt>getAttachedConvTester()->convStatus()<tt> if
	 * <tt>getAttachedConvTester().get()!=NULL &&
	 * attachmentMode()!=ATTACH_TEST_EXCLUDE</tt>.
	 *
	 * The way that that status test returned form
	 * <tt>getAttachedConvTester()->convStatus()<tt> is dealt with is
	 * determined by the value returned by <tt>attachmentMode()</tt>
	 * just before this function is called.  It should be obvious what
	 * the post conditions are base on the value of
	 * <tt>attachmentMode()</tt>.
	 */
	void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);
	///
	void attach( const Teuchos::RefCountPtr<ConvergenceTester<Scalar> > &convTester );
	///
	Teuchos::RefCountPtr<ConvergenceTester<Scalar> > getAttachedConvTester();

	//@}

protected:

	/** @name Protected pure virtual functions to be overridden */
	//@{

	///
	virtual void protectedReset() = 0;

	///
	virtual void protectedConvStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		) = 0;

	//@}

private:

	Teuchos::RefCountPtr<ConvergenceTester<Scalar> > attachedConvTester_;

}; // class AttachConvergenceTesterBase

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_ATTACH_CONVERGENCE_TESTER_BASE_DECL_HPP
