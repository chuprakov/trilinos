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

// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblem.hpp

#ifndef EPETRA_NONLINEAR_PROBLEM_HPP
#define EPETRA_NONLINEAR_PROBLEM_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

namespace Epetra {

///
/** Epetra namespace for a nonlinear problem.
 *
 * ToDo: Finish documentation!
 */
class NonlinearProblem {
public:

  ///
  typedef double Scalar;

  ///
  virtual ~NonlinearProblem() {}

	/// Value for an infinite bound.
	static Scalar infiniteBound();
	
	/** @name Initialization */
	//@{
	
	///
	virtual void initialize( bool testSetup = false ) = 0;
	///
	virtual bool isInitialized() const = 0;
	
	//@}
	
	/** @name Basic information */
	//@{
	
	///
	virtual int Nu() const;
	
	///
	virtual int numResponseFunctions() const;

	//@}

	/** @name VectorSpaces */
	//@{

	///
	virtual Teuchos::RefCountPtr<const Epetra_Map>  map_y() const = 0;
	///
	virtual Teuchos::RefCountPtr<const Epetra_Map>  map_u(int l) const;
	///
	virtual Teuchos::RefCountPtr<const Epetra_Map>  map_c() const = 0;
	///
	virtual Teuchos::RefCountPtr<const Epetra_Map>  map_g() const;

	//@}

	/** @name Bounds */
	//@{

	///
	virtual const Epetra_Vector& yL() const = 0;
	///
	virtual const Epetra_Vector& yU() const = 0;
	///
	virtual const Epetra_Vector& uL(int l) const;
	///
	virtual const Epetra_Vector& uU(int l) const;
	///
	virtual const Epetra_Vector& gL() const;
	///
	virtual const Epetra_Vector& gU() const;

	//@}

	/** @name Initial values (guesses) for state and auxiliary variables */
	//@{

	///
	virtual const Epetra_Vector& y0() const = 0;
	///
	virtual const Epetra_Vector& u0(int l) const;
	//@}

	/** @name Calculation methods */
	//@{

	///
	virtual void calc_c(
		const Epetra_Vector     &y
		,const Epetra_Vector*   u[]
    ,Epetra_Vector          *c
		) const = 0;

	///
	virtual void calc_g(
		const Epetra_Vector     &y
		,const Epetra_Vector*   u[]
    ,Epetra_Vector          *g
		) const;

	//@}

	/** @name Reporting of the final soltuion */
	//@{

	///
	virtual void reportFinalSolution(
		const Epetra_Vector     &y
		,const Epetra_Vector*   u[]
		,bool                   solved
		);

	//@}

}; // class NonlinearProblem

} // namespace Epetra

#endif // EPETRA_NONLINEAR_PROBLEM_HPP

