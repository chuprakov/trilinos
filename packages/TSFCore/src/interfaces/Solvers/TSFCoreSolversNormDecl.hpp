// ///////////////////////////////////////////////////////////////
// TSFCoreSolversNormDecl.hpp

#ifndef TSFCORE_SOLVERS_NORM_DECL_HPP
#define TSFCORE_SOLVERS_NORM_DECL_HPP

#include "TSFCoreSolversTypes.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Abstract interface for computing a norm of a vector.
 *
 * Since the method <tt>norm()</tt> has a good default
 * implementation based on the natural definition of norm,
 * this class can actually be used to instantiate objects.
 *
 * <b>Notes to subclass developers</b>
 *
 * Subclasses should override the multi-vector version
 * <tt>norms()</tt> and leave the default implementation of
 * <tt>norm()</tt> since it calls the multi-vector version.  In any
 * case the vector and multi-vector versions <tt>norm()</tt> and
 * <tt>norms()</tt> must compute the same result.
 */
template<class Scalar>
class Norm {
public:

	///
	/** Compute a norm of a vector.
	 * 
	 * @param  x  [in] The vector for which who's norm is computed.
	 *
	 * @return  Returns a norm of <tt>x</tt>.
	 *
	 * The default implementation calls the mult-vector version <tt>norms()</tt>.
	 */
	virtual Scalar norm(const Vector<Scalar>& x) const;

	///
	/** Compute a norms of each column in a multi-vector.
	 * 
	 * @param  X      [in] The multi-vector for which who's column norms are computed.
	 * @param  norms  [out] Array (length <tt>X.domain()->dim()</tt>) of the norms
	 *                 <tt>norms[j-1] = this->norm(*X.col(j))</tt>, for <tt>j = 1 ... X.domain()->dim()</tt>.
	 *
	 * Postconditions:<ul>
	 * <li><tt>norms[j-1] == this->norm(*X.col(j))</tt>, for <tt>j = 1 ... X.domain()->dim()</tt>
	 * </ul>
	 *
	 * The default implementation computes:
	 * <tt>norms[j-1] = sqrt(scalar_prod[j-1]),</tt>, for <tt>j = 1 ... X.domain()->dim()</tt>.
	 * where <tt>scalar_prod</tt> is a array computed by:
	 * <tt>X.range()->scalarProd(X,X,scalar_prod)</tt>.
	 */
	virtual void norms( const MultiVector<Scalar>& X, Scalar norms[] ) const;

}; // class Norm

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORM_DECL_HPP
