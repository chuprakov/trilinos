// ////////////////////////////////////////////////////////////
// TSFCoreNonlinNP2DSimDecl.hpp

#ifndef TSFCORE_NONLIN_NP_2D_SIM_DECL_HPP
#define TSFCORE_NONLIN_NP_2D_SIM_DECL_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolveIter.hpp"

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
 * is used to solve the linear systems (without a preconditioner).
 */
template<class Scalar>
class NP2DSim : public NonlinearProblemFirstOrder<Scalar> {
public:

	/** @name Constructors / Initializers / accessors */
	//@{

	///
	/** Constructor.
	 *
	 * ToDo: Finish documentation!
	 */
	NP2DSim(
		const Scalar                                                  a           = 2.0
		,const Scalar                                                 b           = 0.0
		,const Scalar                                                 d         = 10.0
		,const Scalar                                                 lin_sol_tol = 1e-12
		,const MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >  &space_y_c  = MemMngPack::null
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
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >  space_y() const;
	///
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >  space_c() const;
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
	MemMngPack::ref_count_ptr< const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > > factory_DcDy() const;
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

	bool                                                                 isInitialized_;
	Scalar                                                               a_;
	Scalar                                                               b_;
	Scalar                                                               d_;
	Scalar                                                               lin_sol_tol_;
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >                space_y_c_;
	MemMngPack::ref_count_ptr<Vector<Scalar> >                           yL_;
	MemMngPack::ref_count_ptr<Vector<Scalar> >                           yU_;
	MemMngPack::ref_count_ptr<Vector<Scalar> >                           y0_;
	MemMngPack::ref_count_ptr<
		const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > >  factory_DcDy_;
	Vector<Scalar>                                                       *c_;
	LinearOpWithSolveIter<Scalar>                                        *DcDy_;

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
