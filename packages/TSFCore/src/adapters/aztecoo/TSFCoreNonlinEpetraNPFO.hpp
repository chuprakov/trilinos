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

// //////////////////////////////////////////////////////
// TSFCoreNonlinEpetraNPFO.hpp

#ifndef TSFCORE_NONLIN_EPETRA_NPFO_HPP
#define TSFCORE_NONLIN_EPETRA_NPFO_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolveAztecOO.hpp"
#include "Epetra_NonlinearProblemFirstOrder.hpp"
#include "Epetra_LinearSystemScaler.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "Ifpack_PrecGenerator.hpp"
#include "TSFCoreLinearOpTester.hpp"
#include "AztecOO.h"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Implements the <tt>TSFCore::Nonlin::NonlinearProblemFirstOrder</tt>
 * interface given a <tt>Epetra::NonlinearProblemFirstOrder</tt> object.
 *
 * Basically, this class wraps Epetra objects in TSFCore/Epetra
 * adapter objects and helps take care of defining algebraic
 * preconditioners with <tt>Ifpack</tt> and solving linear systems
 * with <tt>AztecOO</tt>.
 *
 * ToDo: Finish documentation!
 */
class EpetraNPFO : public NonlinearProblemFirstOrder<double> {
public:

	///
	typedef double Scalar;
	///
	using NonlinearProblem<Scalar>::get_c;
	///
	using NonlinearProblem<Scalar>::get_g;
	///
	using NonlinearProblemFirstOrder<Scalar>::get_DcDy;
	///
	using NonlinearProblemFirstOrder<Scalar>::get_DcDu;
	///
	using NonlinearProblemFirstOrder<Scalar>::get_DgDy;
	///
	using NonlinearProblemFirstOrder<Scalar>::get_DgDu;

	/// Stream that trace to which information will be sent
	STANDARD_NONCONST_COMPOSITION_MEMBERS( std::ostream, trace_out )

	/// Set if the constraints c are to be automatically scaled by inverse row-sums of the initial Jacobian
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, autoScaleStateConstraints )

	/// Set if the state variables y are to be automatically scaled by inverse column-sums of the initial Jacobian
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, autoScaleStateVariables )

  /// Set the maximum number of linear solver iterations
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxLinSolveIter ) 

  /// Set the relative residual tolerance for the linear solver
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, relLinSolveTol )

  /// Determine if preconditioning is used or not
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, usePrec )

  /// Determines if operators are tested after they are formed
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, testOperators )

  /// Determine the file that the initial state solution is read from.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, yGuessFileNameBase )

  /// Determine the file that the final solution for the state is written.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, yFinalFileNameBase )

  ///
  /** Give mutable access to the object used to generate
   * preconditioners.
   *
   * The purpose of this function is to allow clients to change the
   * options that affect how preconditioners are generated.
   */
  Ifpack::PrecGenerator& precGenerator();

  ///
  const Ifpack::PrecGenerator& precGenerator() const;

	///
	/** Give non-const access to the linear system scaler object so that
	 * clients can set options.
	 */
	Epetra::LinearSystemScaler& linearSystemScaler();

	///
	const Epetra::LinearSystemScaler& linearSystemScaler() const;

  ///
  /** Give mutable access to the object used to test the
   * operators for DcDy, DcDu(l).
   *
   * The purpose of this function is to allow clients to change the
   * options that affect how these tests are performed.
   */
  LinearOpTester<Scalar>& linearOpTester();

  ///
  const LinearOpTester<Scalar>& linearOpTester() const;
	
	///
	/** Give mutable access to the AztecOO object that contains the
	 * Aztec options what will be used to solve linear systems.
	 *
	 * Note: Only the Aztec options in the object will be copied off and
	 * nothing else!
	 */
	AztecOO& aztecOO();

	/// Scaling vector for constraints c(y,u)
	STANDARD_COMPOSITION_MEMBERS( EpetraVector, c_scaling )

	/// Scaling vector for state variables y
	STANDARD_COMPOSITION_MEMBERS( EpetraVector, y_scaling )

	/** @name Constructors / Initializers / accessors */
	//@{

  ///
	/** Construct to uninitialized (but with default options set)
	 *
	 * Note, these defaults where taken from
	 * NOX::Epetra::Group::applyJacobianInverse(...) on 2004/01/19.
	 */
  EpetraNPFO(
		const bool           autoScaleStateConstraints = false
		,const bool          autoScaleStateVariables   = false
	 	,const int           maxLinSolveIter      = 400
		,const double        relLinSolveTol       = 1e-6
		,const bool          usePrec              = true
		,const bool          testOperators        = false
		,const std::string   &yGuessFileNameBase  = ""
		,const std::string   &yFinalFileNameBase  = ""
		);

  ///
	/** Initialize this nonlinear problem.
	 *
	 * ToDo: Finish documentation!
	 */
  void initialize(
    const Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>   &epetra_np
    );

	//@}

	/** @name Overridden from NonlinearProblem */
	//@{

	///
	void initialize( bool testSetup );
	///
	bool isInitialized() const;
	///
	int Nu() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_y() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_u(int l) const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_c() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_g() const;
	///
	const Vector<Scalar>& yL() const;
	///
	const Vector<Scalar>& yU() const;
	///
	const Vector<Scalar>& uL(int l) const;
	///
	const Vector<Scalar>& uU(int l) const;
	///
	const Vector<Scalar>& gL() const;
	///
	const Vector<Scalar>& gU() const;
	///
	const Vector<Scalar>& y0() const;
	///
	const Vector<Scalar>& u0(int l) const;
	///
	void set_c(Vector<Scalar>* c);
	///
	Vector<Scalar>* get_c();
	///
	void set_g(Vector<Scalar>* g);
	///
	Vector<Scalar>* get_g();
	///
	void unsetQuantities();
	///
	void calc_c(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_g(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void reportFinalSolution(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    solved
		);
	
	//@}

	/** @name Overridden from NonlinearProblemFirstOrder */
	//@{

	///
	Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOpWithSolve<Scalar> > > factory_DcDy() const;
	///
	Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOp<Scalar > > > factory_DcDu(int l) const;
	///
	ETransp opDcDy() const;
	///
	ETransp opDcDu(int l) const;
	///
	void set_DcDy(LinearOpWithSolve<Scalar>* DcDy);
	///
	LinearOpWithSolve<Scalar>* get_DcDy();
	///
	void set_DcDu(int l, LinearOp<Scalar>* DcDu_l);
	///
	LinearOp<Scalar>* get_DcDu(int l);
	///
	void set_DgDy(MultiVector<Scalar>* DgDy);
	///
	MultiVector<Scalar>* get_DgDy();
	///
	void set_DgDu(int l, MultiVector<Scalar>* DgDu_l);
	///
	MultiVector<Scalar>* get_DgDu(int l);
	///
	void calc_DcDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_DcDu(
		int                      l
		,const Vector<Scalar>    &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_DgDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_DgDu(
		int                      l
		,const Vector<Scalar>    &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;

	//@}

private:

  // ///////////////////////////////////////
  // Private types

public: // Intel C++ 6.0 requires this?

  ///
  class DcDu_Allocator {
  public:
    
    /** @name Constructors */
    //@{
    
    ///
    DcDu_Allocator(
      const bool  useEO
      );
      
    //@}
      
    /** @name AbstractFactoryStd compliant interface */
    //@{
    
    ///
    typedef Teuchos::RefCountPtr<TSFCore::LinearOp<double> >  ptr_t;
    ///
    const ptr_t allocate() const;

    //@}
    
  private:
    
    DcDu_Allocator();  // Not derfined and not to be called
    
    const bool  useEO_;
    
  }; // class DcDu_Allocator

private:

	// //////////////////////////////////////
	// Private data members

	bool isInitialized_;

	int numProc_;
	int procRank_;
	
	mutable AztecOO  aztecOO_;

  Ifpack::PrecGenerator        precGenerator_;
	Epetra::LinearSystemScaler   linearSystemScaler_;
	LinearOpTester<Scalar>       linearOpTester_;

  Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>        epetra_np_;

	Teuchos::RefCountPtr<const EpetraVectorSpace >                  space_y_;
	std::vector<Teuchos::RefCountPtr<const EpetraVectorSpace > >    space_u_;
	Teuchos::RefCountPtr<const EpetraVectorSpace >                  space_c_;
	Teuchos::RefCountPtr<const EpetraVectorSpace >                  space_g_;

	Teuchos::RefCountPtr<EpetraVector>                           yL_;
	Teuchos::RefCountPtr<EpetraVector>                           yU_;
	Teuchos::RefCountPtr<EpetraVector>                           y0_;
	std::vector<Teuchos::RefCountPtr<EpetraVector> >             uL_;
 	std::vector<Teuchos::RefCountPtr<EpetraVector> >             uU_;
	std::vector<Teuchos::RefCountPtr<EpetraVector> >             u0_;
	Teuchos::RefCountPtr<EpetraVector>                           gL_;
 	Teuchos::RefCountPtr<EpetraVector>                           gU_;

	Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolve<Scalar> > >       factory_DcDy_;
  std::vector<Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOp<Scalar> > > >  factory_DcDu_;

	Teuchos::RefCountPtr<EpetraVector>        y_scaling_inv_;

	mutable std::vector<const Epetra_Vector*>  u_unscaled_;
	Teuchos::RefCountPtr<Epetra_Vector>        y_unscaled_;

  mutable bool c_updated_, g_updated_, DcDy_updated_, DgDy_updated_;
  mutable std::vector<bool> DcDu_updated_, DgDu_updated_;

  EpetraVector                                    *c_;
	EpetraVector                                    *g_;
	LinearOpWithSolveAztecOO                        *DcDy_;
  std::vector<EpetraLinearOp*>                    DcDu_op_;
  std::vector<EpetraMultiVector*>                 DcDu_mv_;
  EpetraMultiVector                               *DgDy_;
  std::vector<EpetraMultiVector*>                 DgDu_;

  mutable std::vector<Teuchos::RefCountPtr<Epetra_Operator> >     epetra_DcDu_op_;
  mutable std::vector<Teuchos::RefCountPtr<Epetra_MultiVector> >  epetra_DcDu_mv_;
  mutable std::vector<Epetra::EpetraOp_or_EpetraMV>               epetra_DcDu_args_;

  mutable std::vector<Teuchos::RefCountPtr<Epetra_MultiVector> >  epetra_DgDu_;
  mutable std::vector<Epetra_MultiVector*>                        epetra_DgDu_args_;

	// //////////////////////////////////////
	// Private member functions

	//
	void read_y_guess( EpetraVector *y );
	//
	void write_y_final( const Epetra_Vector &y );
  //
  const Epetra_Vector& set_y( const Vector<Scalar> &y_scaled ) const;
  //
	const Epetra_Vector** set_u( const Vector<Scalar>* u[], bool newPoint ) const;
	//
	void computeScaling();
	//
	void scale_y( const Epetra_MultiVector &y_unscaled, Epetra_MultiVector *y_scaled ) const;
	//
	void unscale_y( const Epetra_MultiVector &y_scaled, Epetra_MultiVector *y_unscaled ) const;
	//
	void scale_c( const Epetra_MultiVector &c_unscaled, Epetra_MultiVector *c_scaled ) const;
  //
	void updateNewPoint( bool newPoint ) const;
  //
  void calc_Dc(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
    ,bool                    computeGradients
    ) const;
	//
  void calc_Dg(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
    ,bool                    computeGradients
    ) const;

	// Not defined and not to be called
	EpetraNPFO(const EpetraNPFO&);
	EpetraNPFO& operator=(const EpetraNPFO&);

}; // class EpetraNPFO

// ///////////////////////////////////
// Inline members

// public

inline
Ifpack::PrecGenerator& EpetraNPFO::precGenerator()
{
  return precGenerator_;
}

inline
const Ifpack::PrecGenerator& EpetraNPFO::precGenerator() const
{
  return precGenerator_;
}

inline
Epetra::LinearSystemScaler&
EpetraNPFO::linearSystemScaler()
{
	return linearSystemScaler_;
}

inline
const Epetra::LinearSystemScaler&
EpetraNPFO::linearSystemScaler() const
{
	return linearSystemScaler_;
}

inline
LinearOpTester<EpetraNPFO::Scalar>& EpetraNPFO::linearOpTester()
{
	return linearOpTester_;
}

inline
const LinearOpTester<EpetraNPFO::Scalar>& EpetraNPFO::linearOpTester() const
{
	return linearOpTester_;
}

inline
AztecOO& EpetraNPFO::aztecOO()
{
	return aztecOO_;
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_EPETRA_NPFO_HPP
