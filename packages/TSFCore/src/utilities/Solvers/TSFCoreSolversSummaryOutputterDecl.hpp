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
// TSFCoreSolversSummaryOutputterDecl.hpp

#ifndef TSFCORE_SOLVERS_SUMMARY_OUTPUTTER_DECL_HPP
#define TSFCORE_SOLVERS_SUMMARY_OUTPUTTER_DECL_HPP

#include "TSFCoreSolversAttachConvergenceTesterBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "StandardCompositionMacros.hpp"

namespace TSFCore {
namespace Solvers {

///
/** A convergence tester subclass that print iteration summary
 * information only.
 *
 * This class will always just applies and returns the convergence
 * test of any attached object (if a client calls <tt>attach()</tt>).
 */
template<class Scalar>
class SummaryOutputter : public AttachConvergenceTesterBase<Scalar> {
public:

	/// Set the stream that output will be sent to
	STANDARD_COMPOSITION_MEMBERS( std::ostream, out )

	/// Set the leading string that will be printed at the beginning of each new line of output.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, leadingOutputStr )

	/** @name Constructors / initializers */
	//@{
	
	///
	/** Construct to uninitialized.
	 */
	SummaryOutputter();

	///
	/** Construct with output stream and leading string optionally set.
	 *
	 * @param  out    [in] Smart pointer to output stream.
	 * @param  leadingOutputStr
	 *                [in] String that is printed at beginning of each new line.
	 *
	 * Preconditions:<ul>
	 * <li><tt></tt>
	 * </ul>
	 */
	SummaryOutputter(
		const out_ptr_t       &out
		,const std::string    &leadingOutputStr
		);

	//@}

	/** @name Overridden from ConvergenceTester */
	//@{

	///
	void reset();
	///
	void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);

	//@}

protected:

	/** @name Overridden from AttachedConvergenceTesterBase */
	//@{

	/// Never called
	void protectedReset();
	/// Never called
	void protectedConvStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);

	//@}

private:

	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  ScalarMagnitude;

	bool resetCalled_;
	std::valarray<ScalarMagnitude> norms_;

}; // class SummaryOutputter

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_SUMMARY_OUTPUTTER_DECL_HPP
