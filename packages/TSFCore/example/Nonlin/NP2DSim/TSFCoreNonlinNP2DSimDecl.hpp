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

// ////////////////////////////////////////////////////////////
// TSFCoreNonlinNP2DSimDecl.hpp

#ifndef TSFCORE_NONLIN_NP_2D_SIM_DECL_HPP
#define TSFCORE_NONLIN_NP_2D_SIM_DECL_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolveIter.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Concrete implementation of <tt>NonlinearProblemFirstOrder</tt> for
 * a simple set of 2-D equations.
 *
 * The simulation equations represented by this class are
 \verbatim

    c(1) =       y(1)   + y(2)^2 - a   = 0
    c(2) = d * ( y(1)^2 - y(2)   - b ) = 0
 \endverbatim
 * which only has one unique solution (1.0,1.0) if
 * <tt>a = 2.0</tt> and <tt>b = 0.0</tt> (which is the
 * default).
 *
 * The Jacobian <tt>DcDy</tt> of the above equations is
 \verbatim

    DcDy = [  1.0       2*y(2)  ]
           [  2*d*y(1)  -d      ]
 \endverbatim
 * This Jacobian is nonsingular for every point except (1/(2*d),-0.5)
 * and (-1/(2*d),+0.5).
 *
 * The vector space <tt>SerialVectorSpace</tt> is used as the default
 * for <tt>space_y()</tt> and <tt>space_c()</tt> and a
 * <tt>MultiVector</tt> is used for the basic <tt>LinearOp</tt> object
 * for the Jacobian <tt>DcDy</tt>.  The class
 * <tt>LinearOpWithSolveIter</tt> is used for the nonsingualar
 * Jacobian object <tt>DcDy</tt> and the subclass <tt>BiCGSolver</tt>
 * is used to solve the linear systems.
 */
template<class Scalar>
class NP2DSim : public NonlinearProblemFirstOrder<Scalar> {
public:

  ///
  //using NonlinearProblem<Scalar>::get_c;
  ///
  using NonlinearProblemFirstOrder<Scalar>::get_DcDy;

	/** @name Public types and options */
	//@{

	///
	enum ELinearSolverType { LINEAR_SOLVER_BICG, LINEAR_SOLVER_GMRES };
	
	/// Set the type of linear solver used.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ELinearSolverType, linearSolverType )	

  /// Determine if we will dump to a file or not.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dumpToFile )	

	//@}

	/** @name Constructors / Initializers / accessors */
	//@{

	///
	/** Constructor.
	 *
	 * ToDo: Finish documentation!
	 */
	NP2DSim(
		const Scalar                                             a                  = 2.0
		,const Scalar                                            b                  = 0.0
		,const Scalar                                            d                  = 10.0
		,const Scalar                                            lin_sol_tol        = 1e-12
		,const ELinearSolverType                                 linearSolverType   = LINEAR_SOLVER_GMRES  // Test GMRES
		,const bool                                              dumpToFile         = false
		,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >  &space_y_c         = Teuchos::null
		);
	///
	void set_a( const Scalar a );
	///
	Scalar get_a() const;
	///
	void set_b( const Scalar b );
	///
	Scalar get_b() const;
	///
	void set_d( const Scalar d );
	///
	Scalar get_d() const;
	///
	void set_y0( const Scalar y01, const Scalar y02 );

	//@}

	/** @name Overridden from NonlinearProblem */
	//@{

	///
	void initialize( bool testSetup );
	///
	bool isInitialized() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_y() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_c() const;
	///
	const Vector<Scalar>& yL() const;
	///
	const Vector<Scalar>& yU() const;
	///
	const Vector<Scalar>& y0() const;
	///
	void set_c(Vector<Scalar>* c);
	///
	Vector<Scalar>* get_c();
	///
	void unsetQuantities();
	///
	void calc_c(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;

	//@}

	/** @name Overridden from NonlinearProblemFirstOrder */
	//@{
	
	///
	Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOpWithSolve<Scalar> > > factory_DcDy() const;
	///
	ETransp opDcDy() const;
	///
	void set_DcDy(LinearOpWithSolve<Scalar>* DcDy);
	///
	LinearOpWithSolve<Scalar>* get_DcDy();
	///
	void calc_DcDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;

	//@}

private:

	bool                                                            isInitialized_;
	Scalar                                                          a_;
	Scalar                                                          b_;
	Scalar                                                          d_;
	Scalar                                                          lin_sol_tol_;
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >                space_y_c_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           yL_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           yU_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           y0_;
	Teuchos::RefCountPtr<
		const Teuchos::AbstractFactory<LinearOpWithSolve<Scalar> > >  factory_DcDy_;
	Vector<Scalar>                                                  *c_;
	LinearOpWithSolveIter<Scalar>                                   *DcDy_;

  mutable Teuchos::RefCountPtr<std::ostream>                      linear_solver_out_;

	// Not defined and not to be called
	NP2DSim(const NP2DSim<Scalar>&);
	NP2DSim<Scalar>& operator=(const NP2DSim<Scalar>&);

}; // class NP2DSim

// //////////////////////////
// Inline members

template<class Scalar>
inline
void NP2DSim<Scalar>::set_a( const Scalar a )
{
	a_ = a;
}

template<class Scalar>
inline
Scalar NP2DSim<Scalar>::get_a() const
{
	return a_;
}

template<class Scalar>
inline
void NP2DSim<Scalar>::set_b( const Scalar b )
{
	 b_ = b;
}

template<class Scalar>
inline
Scalar NP2DSim<Scalar>::get_b() const
{
	return b_;
}

template<class Scalar>
inline
void NP2DSim<Scalar>::set_d( const Scalar d )
{
	 d_ = d;
}

template<class Scalar>
inline
Scalar NP2DSim<Scalar>::get_d() const
{
	return d_;
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NP_2D_SIM_DECL_HPP
