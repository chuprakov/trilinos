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
// TSFCoreNonlinNP4DOptDecl.hpp

#ifndef TSFCORE_NONLIN_NP_4D_OPT_DECL_HPP
#define TSFCORE_NONLIN_NP_4D_OPT_DECL_HPP

#include "TSFCoreNonlinNP2DSim.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Concrete implementation of <tt>NonlinearProblemFirstOrder</tt> for
 * a simple set of 2-D equations with 2-D design variables.
 *
 * This subclass shows how it is possible to build on an existing
 * <tt>NonlinearProblemFirstOrder</tt> class (which was developed for
 * simulation only) and to augment it with with auxiliary variables
 * and response functions.  In this case, the equations modeled in
 * <tt>NP2DSim</tt> are augmented with two auxiliary variables and
 * four response functions.
 *
 * This subclass represents the following state equations and response
 * functions:
 \verbatim

    c(1)  =    y(1)   + y(2)^2 - u(1)  = 0
    c(2)  = d*(y(1)^2 - y(2)   - u(2)) = 0
    -10  <=    g(1) = y(1) - yt1       <= +10
    -10  <=    g(2) = y(2) - yt2       <= +10
    -10  <=    g(3) = u(1) - ut1       <= +10
    -10  <=    g(4) = u(2) - ut2       <= +10
 \endverbatim
 * Above, <tt>u(1)</tt> and <tt>u(2)</tt> are substituted into
 * <tt>a</tt> and <tt>b</tt> in the model equations in
 * <tt>NP2DSim</tt>.  Note that the above set of state equations and
 * response functions can represent a square 4-D set of equations
 * (with <tt>g(y,u)(1,2) == 0</tt>) or a composite objective function
 * can formed from <tt>g(y,u)</tt> (i.e. <tt>f(y,u) = sum(
 * 0.5*(g(k)^2), for k = 1..4 )</tt>).  This subclass does not fix how
 * the auxiliary response functions are treated.
 *
 * The initial guess <tt>y0</tt> is determined by the underlying
 * <tt>NP2DSim</tt> object and the initial guess <tt>u0</tt> is
 * <tt>[ut1,ut2]</tt>.
 *
 * ToDo: Finish documentation
 */
template<class Scalar>
class NP4DOpt : public NonlinearProblemFirstOrder<Scalar> {
public:

  ///
  //using NonlinearProblem<Scalar>::get_c;
  ///
  //using NonlinearProblem<Scalar>::get_g;
  ///
  using NonlinearProblemFirstOrder<Scalar>::get_DcDy;
	///
  using NonlinearProblemFirstOrder<Scalar>::get_DcDu;
	///
  using NonlinearProblemFirstOrder<Scalar>::get_DgDy;
	///
  using NonlinearProblemFirstOrder<Scalar>::get_DgDu;

	/** @name Constructors / Initializers / accessors */
	//@{

	///
	/** Constructor.
	 *
	 * ToDo: Finish documentation!
	 */
	NP4DOpt(
		const Scalar                                             yt1                = 1.0
		,const Scalar                                            yt2                = 1.0
		,const Scalar                                            ut1                = 2.0
		,const Scalar                                            ut2                = 0.0
		,const Scalar                                            d                  = 10.0
		,const Scalar                                            lin_sol_tol        = 1e-12
		,const typename NP2DSim<Scalar>::ELinearSolverType       linearSolverType   = NP2DSim<Scalar>::LINEAR_SOLVER_BICG // Test BiCG
		,const bool                                              dumpToFile         = false
		,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >  &space_y_c         = Teuchos::null
		);
	///
	void set_ut1( const Scalar ut1 );
	///
	Scalar get_ut1() const;
	///
	void set_ut2( const Scalar ut2 );
	///
	Scalar get_ut2() const;
	///
	void set_yt1( const Scalar yt1 );
	///
	Scalar get_yt1() const;
	///
	void set_yt2( const Scalar yt2);
	///
	Scalar get_yt2() const;
	///
	void set_y0( const Scalar y01, const Scalar y02 );
	///
	void set_u0( const Scalar u01, const Scalar u02 );

	//@}

	/** @name Overridden from NonlinearProblem */
	//@{

	///
	void initialize( bool testSetup );
	///
	bool isInitialized() const;
	/// Returns 1
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

	//@}

	/** @name Overridden from NonlinearProblemFirstOrder */
	//@{

	///
	Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > > factory_DcDy() const;
	///
	Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOp<Scalar > > > factory_DcDu(int l) const;
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

	enum { NUM_RESPONSE_FUNCTIONS = 4 };

	// //////////////////////////////////////
	// Private data members

	NP2DSim<Scalar>                                                      np2dsim_;
	Scalar                                                               yt1_;
	Scalar                                                               yt2_;
	Scalar                                                               ut1_;
	Scalar                                                               ut2_;
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >                space_g_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           uL_;
 	Teuchos::RefCountPtr<Vector<Scalar> >                           uU_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           gL_;
 	Teuchos::RefCountPtr<Vector<Scalar> >                           gU_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           u0_;
	Teuchos::RefCountPtr<
		const MemMngPack::AbstractFactory<LinearOp<Scalar> > >           factory_DcDu_;

	Vector<Scalar>                                                       *g_;
    MultiVector<Scalar>                                                  *DcDu_;
    MultiVector<Scalar>                                                  *DgDy_;
    MultiVector<Scalar>                                                  *DgDu_;

	// //////////////////////////////////////
	// Private member functions

	void set_u(const Vector<Scalar>* u[], bool newPoint) const;

	// Not defined and not to be called
	NP4DOpt(const NP4DOpt<Scalar>&);
	NP4DOpt<Scalar>& operator=(const NP4DOpt<Scalar>&);

}; // class NP4DOpt

// //////////////////////////
// Inline members

template<class Scalar>
inline
void NP4DOpt<Scalar>::set_ut1( const Scalar ut1 )
{
	ut1_ = ut1;
}

template<class Scalar>
inline
Scalar NP4DOpt<Scalar>::get_ut1() const
{
	return ut1_;
}

template<class Scalar>
inline
void NP4DOpt<Scalar>::set_ut2( const Scalar ut2 )
{
	 ut2_ = ut2;
}

template<class Scalar>
inline
Scalar NP4DOpt<Scalar>::get_ut2() const
{
	return ut2_;
}

template<class Scalar>
inline
void NP4DOpt<Scalar>::set_yt1( const Scalar yt1 )
{
	 yt1_ = yt1;
}

template<class Scalar>
inline
Scalar NP4DOpt<Scalar>::get_yt1() const
{
	return yt1_;
}


template<class Scalar>
inline
void NP4DOpt<Scalar>::set_yt2( const Scalar yt2 )
{
	 yt2_ = yt2;
}

template<class Scalar>
inline
Scalar NP4DOpt<Scalar>::get_yt2() const
{
	return yt2_;
}

template<class Scalar>
inline
void NP4DOpt<Scalar>::set_y0( const Scalar y01, const Scalar y02 )
{
	np2dsim_.set_y0(y01,y02);
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NP_4D_OPT_DECL_HPP
