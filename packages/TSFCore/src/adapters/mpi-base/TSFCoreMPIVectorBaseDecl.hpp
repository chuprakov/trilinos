// /////////////////////////////////////////////////////////////////
// TSFCoreMPIVectorBaseDecl.hpp

#ifndef TSFCORE_MPI_VECTOR_BASE_DECL_HPP
#define TSFCORE_MPI_VECTOR_BASE_DECL_HPP

#include "TSFCoreVector.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

///
/** Interface base subclass for  MPI-based vectors.
 *
 * By inheriting from this base class, vector implementations allow
 * their vector objects to be seamlessly combined with other MPI-based
 * vector objects (of different concrete type) in <tt>applyOp()</tt>.
 * A big part of this protocal is that every vector object can expose
 * an <tt>MPIVectorSpaceBase<></tt> object through the virtual
 * <tt>mpiSpace()</tt> method.
 *
 * This node subclass contains an implementation of
 * <tt>applyOp()</tt> that relies on implementations of the methods
 * (<tt>const</tt>) <tt>getSubVector()</tt>, <tt>freeSubVector()</tt>,
 * (non-<tt>const</tt>) <tt>getSubVector()</tt> and
 * <tt>commitSubVector()</tt>.  In essense, this implemenation will
 * only call the <tt>getSubVector()</tt> methods using a range of
 * (global) indexes for elements that exist on the local processor.
 * As long as the number of local elements on each processor is fairly
 * large, the virtual function call overhead will be minimal and this
 * will result in a near optimal implementation.
 *
 * The only methods that a subclass must override in order to
 * implement a fully-functional vector subclass of which who's vector
 * objects can be seemlessly used with any other MPI-based vector
 * objects are the methods (<tt>const</tt>) <tt>getSubVector()</tt>,
 * <tt>freeSubVector()</tt>, (non-<tt>const</tt>)
 * <tt>getSubVector()</tt> and <tt>commitSubVector()</tt> and
 * <tt>mpiSpace()</tt>.
 *
 * If the <tt>getSubVector()</tt> methods are ever called with index
 * ranges outside of those of the local processor, then the default
 * implementations of all of the methods (<tt>const</tt>)
 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt>,
 * (non-<tt>const</tt>) <tt>getSubVector()</tt> and
 * <tt>commitSubVector()</tt> should be called in instead (see the
 * vector subclass <tt>EpetraVector</tt> for an example).
 * Alternatively, a subclass could provide more specialized
 * implemenations of these methods (for more efficient gather/scatter
 * operations) if desired but this should not be needed for most use
 * cases.
 *
 * It is interesting to note that is case the explicit subvector
 * access methods call on its default implementation defined in
 * <tt>Vector</tt> (which calls on <tt>applyOp()</tt>), the operator
 * will be properly applied since the version of <tt>applyOp()</tt>
 * implemented in this class will only request local vector data and
 * hence there will only be two levels of recussion for any call to an
 * explicit subvector access method.  This is a truly elegat result.
 *
 * As described in the documentation for <tt>MPIVectorSpaceBase</tt>
 * it is possible that at runtime that it may be discovered that
 * mapping of vector data to processors does not fall under this
 * design in which case the method <tt>applyOp()</tt> should be
 * overridden to handle this which will of course remove the
 * possibility of interoperability with other MPI-based vector
 * objects.
 */
template<class Scalar>
class MPIVectorBase : virtual public Vector<Scalar> {
public:

	///
	MPIVectorBase();

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI-based vector space object for <tt>*this</tt> vector.
	 */
	virtual MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

	//@}

	/** @name Overridden from Vector */
	//@{

	/// Calls <tt>mpiSpace()</tt>
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > space() const;

	///
	/** Implements the <tt>%applyOp()</tt> method through the methods
	 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt> and
	 * <tt>commitSubVector()</tt> as described above.
	 *
	 * Note that if this method is entered again before a call has
	 * been completed, then this is an indication that the methods
	 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt> and/or
	 * <tt>commitSubVector()</tt> have not been overridden properly
	 * and this method will then throw an exception.
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

}; // end class MPIVectorBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_BASE_DECL_HPP
