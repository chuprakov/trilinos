// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreOpBaseDecl.hpp

#ifndef TSFCORE_OP_BASE_DECL_HPP
#define TSFCORE_OP_BASE_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

///
/** Base class for all operators.
 *
 * It is not expected that clients will manipulate objects directly through this
 * interface.  This interface is just ment to provide common declarations
 * the methods <tt>domain()</tt>, <tt>range()</tt> and <tt>opSupported()</tt>.
 */
template<class Scalar>
class OpBase {
public:

	///
	virtual ~OpBase();

	/** @name Pure virtual methods (must be overridden by subclass) */

	///
	/** Domain space for <tt>this</tt> operator.
	 *
	 * Note that a return value of <tt>return.get()==NULL</tt> is a flag that <tt>*this</tt>
	 * is not fully initialized.
	 *
	 * If <tt>return.get()!=NULL</tt>, it is required that the object
	 * referenced by <tt>*return.get()</tt> must have lifetime that
	 * extends past the lifetime of the returned smart pointer object.
	 * However, the object referenced by <tt>*return.get()</tt> my
	 * change if <tt>*this</tt> modified so this reference should not
	 * be maintained for too long.
	 */
	virtual MemMngPack::ref_count_ptr< const VectorSpace<Scalar> > domain() const = 0;

	///
	/** Range space for <tt>this</tt> operator.
	 *
	 * Note that a return value of <tt>return.get()==NULL</tt> is a flag that <tt>*this</tt>
	 * is not fully initialized.
	 *
	 * If <tt>return.get()!=NULL</tt>, it is required that the object
	 * referenced by <tt>*return.get()</tt> must have lifetime that
	 * extends past the lifetime of the returned smart pointer object.
	 * However, the object referenced by <tt>*return.get()</tt> my
	 * change if <tt>*this</tt> modified so this reference should not
	 * be maintained for too long.
	 */
	virtual MemMngPack::ref_count_ptr< const VectorSpace<Scalar> > range() const = 0;

	//@}

	/** @name Virtual functions with default implemenations */
	//@{

	///
	/** Return if the <tt>M_trans</tt> operation is supported.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>true</tt>.
	 *
	 * Note that an operator must support at least one of the values
	 * of <tt>ETrans</tt> (i.e. the transposed or the nontranspoed
	 * operations must be supported, both can not be unsupported)
	 */
	virtual bool opSupported(ETransp M_trans) const;

	//@}

#ifdef DOXYGEN_COMPILE
	const VectorSpace<Scalar>*  domain; // doxygen only!
	const VectorSpace<Scalar>*  range;  // doxygen only!
#endif

};	// end class OpBase

}	// end namespace TSFCore

#endif	// TSFCORE_OP_BASE_DECL_HPP
