// /////////////////////////////////////////////////////////////////
// SerialVectorBaseDecl.hpp

#ifndef TSFCORE_VECTOR_SERIAL_BASE_DECL_HPP
#define TSFCORE_VECTOR_SERIAL_BASE_DECL_HPP

#include "TSFCoreVectorDecl.hpp"

namespace TSFCore {

///
/** Node subclass of serial vectors.
 *
 * This node subclass contains the an implementation of
 * <tt>applyOp()</tt> that relies on implementations of the methods
 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt> and
 * <tt>commitSubVector()</tt>.  A concrete subclass must
 * implement these methods without relying on <tt>applyOp()</tt>
 * (see the concrete subclass <tt>SerialVector</tt>).
 */
template<class Scalar>
class SerialVectorBase : virtual public Vector<Scalar> {
public:

	///
	SerialVectorBase();

	/** @name Overridden from Vector */
	//@{

	///
	/** Implements this method through the methods
	 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt> and
	 * <tt>commitSubVector()</tt>.
	 *
	 * Note that if this method is entered again before a call has
	 * been completed, then this is an indication that the methods
	 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt> and/or
	 * <tt>commitSubVector()</tt> have not been overridden properly.
	 */
	void applyOp(
		const RTOpPack::RTOpT<Scalar>   &op
		,const size_t                   num_vecs
		,const Vector<Scalar>*          vecs[]
		,const size_t                   num_targ_vecs
		,Vector<Scalar>*                targ_vecs[]
		,RTOp_ReductTarget              reduct_obj
		,const Index                    first_ele
		,const Index                    sub_dim
		,const Index                    global_offset
		) const;

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	mutable bool in_applyOp_;

}; // end class SerialVectorBase

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_BASE_DECL_HPP
