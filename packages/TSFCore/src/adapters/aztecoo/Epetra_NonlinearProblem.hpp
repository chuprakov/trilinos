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

}; // class NonlinearProblem

} // namespace Epetra

#endif // EPETRA_NONLINEAR_PROBLEM_HPP

