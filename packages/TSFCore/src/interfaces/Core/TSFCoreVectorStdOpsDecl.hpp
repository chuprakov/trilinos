// ////////////////////////////////////////////////////////////////////
// VectorStdOpsDecl.hpp

#ifndef TSFCORE_VECTOR_STD_OPS_DECL_HPP
#define TSFCORE_VECTOR_STD_OPS_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

/** \defgroup TSFCore_VectorStdOps_grp Collection of standard vector operations.
 */
//@{

///
/** <tt>result = sum( v(i), i = 1...v.space()->dim() )</tt>
 */
template<class Scalar>
Scalar sum( const Vector<Scalar>& v );

///
/** <tt>result = ||v||1</tt>
 */
template<class Scalar>
Scalar norm_1( const Vector<Scalar>& v );

///
/** <tt>result = ||v||2</tt>
 */
template<class Scalar>
Scalar norm_2( const Vector<Scalar>& v );

///
/** <tt>result = ||v||inf</tt>
 */
template<class Scalar>
Scalar norm_inf( const Vector<Scalar>& v_rhs );

///
/** <tt>result = x'*y</tt>
 */
template<class Scalar>
Scalar dot( const Vector<Scalar>& x, const Vector<Scalar>& y );

///
/** <tt>result = v(i)</tt>
 */
template<class Scalar>
Scalar get_ele( const Vector<Scalar>& v, Index i );

///
/** <tt>v(i) = alpha</tt>
 */
template<class Scalar>
void set_ele( Index i, Scalar alpha, Vector<Scalar>* v );

///
/** <tt>y = alpha</tt>
 */
template<class Scalar>
void assign( Vector<Scalar>* y, const Scalar& alpha );

///
/** <tt>y = x</tt>
 */
template<class Scalar>
void assign( Vector<Scalar>* y, const Vector<Scalar>& x );

///
/** <tt>y += alpha</tt>
 */
template<class Scalar>
void Vp_S( Vector<Scalar>* y, const Scalar& alpha );

///
/** <tt>y *= alpha</tt>
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>y = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 */
template<class Scalar>
void Vt_S( Vector<Scalar>* y, const Scalar& alpha );

///
/** <tt>y = alpha * x + y</tt>
 */
template<class Scalar>
void Vp_StV( Vector<Scalar>* y, const Scalar& alpha, const Vector<Scalar>& x );

///
/** <tt>y(i) += alpha * x(i) * v(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void ele_wise_prod( const Scalar& alpha, const Vector<Scalar>& x, const Vector<Scalar>& v, Vector<Scalar>* y );

///
/** <tt>y(i) = alpha * x(i) / v(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void ele_wise_divide( const Scalar& alpha, const Vector<Scalar>& x, const Vector<Scalar>& v, Vector<Scalar>* y );

///
/** Seed the random number generator used in <tt>random_vector</tt>
 */
template<class Scalar>
void seed_randomize( unsigned int );

///
/** Generate a random vector with elements uniformly distrubuted
 * elements.
 * 
 * The elements <tt>v->getEle(i)</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using <tt>seed_randomize()</tt>
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, Vector<Scalar>* v );

//@}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_STD_OPS_DECL_HPP
