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
// NormedConvergenceTesterDecl.hpp

#ifndef TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
#define TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP

#include "TSFCoreSolversConvergenceTester.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Convergence test based on relative errors.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NormedConvergenceTester : public ConvergenceTester<Scalar> {
public:

	/** @name Constructors / initializers */
	//@{
	
	/// Calls <tt>initialize()</tt>
	NormedConvergenceTester(
		const Scalar                                           tol             = 1e-12
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	/// Calls <tt>initialize()</tt>
	NormedConvergenceTester(
		const Index                                            totalNumSystems
		,const Scalar                                          tols[]
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	///
	/** Initialize with one set of tolerances for multiple linear systems
	 */
	void initialize(
		const Scalar                                           tol             = 1e-12
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	///
	/** Initialize with different tolerances for multiple linear systems.
	 */
	void initialize(
		const Index                                            totalNumSystems
		,const Scalar                                          tols[]
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	/// Return the minimum (over all iterations) maximum (over all right-hand-sides) tolerance seen
	Scalar minMaxError() const;

	//@}

	/** @name Overridden from ConvergenceTester */
	//@{

	///
	Teuchos::RefCountPtr<const Norm<Scalar> > norm() const;
	///
	void reset();
	///
	void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);

	//@}

private:

	Teuchos::RefCountPtr<const Norm<Scalar> >   norm_;
	std::vector<Scalar>                              tols_;           // if size()==1 then works for multiple systems
	std::vector<Index>                               activeSystems_;  // cache
	std::vector<Scalar>                              norms_;          // cache

	Scalar                                           minMaxErr_;

}; // class NormedConvergenceTester

// ///////////////////////////
// Inline members

template<class Scalar>
inline
Scalar NormedConvergenceTester<Scalar>::minMaxError() const
{
	return minMaxErr_;
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
