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
 * vector objects (of different concrete types) in <tt>applyOp()</tt>.
 * A big part of this protocol is that every vector object can expose
 * an <tt>MPIVectorSpaceBase</tt> object through the virtual function
 * <tt>mpiSpace()</tt>.
 *
 * This node subclass contains an implementation of <tt>applyOp()</tt>
 * that relies on implementations of the methods (<tt>const</tt>)
 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt>,
 * (non-<tt>const</tt>) <tt>getSubVector()</tt> and
 * <tt>commitSubVector()</tt> (which all have default implementations
 * in this subclass).  In essense, this implemenation will only call
 * the <tt>getSubVector()</tt> methods using a range of (global)
 * indexes for elements that exist on the local processor.  As long as
 * the number of local elements on each processor is fairly large, the
 * virtual function call overhead will be minimal and this will result
 * in a near optimal implementation.
 *
 * The only methods that a subclass must override in order to
 * implement a fully-functional vector subclass of which who's vector
 * objects can be seemlessly used with any other MPI-based vector
 * object are the methods <tt>getLocalData()</tt> and
 * <tt>mpiSpace()</tt>.  All of the other methods have good default
 * implementations.  Providing an implementation for <tt>mpiSpace</tt>
 * of course means having the define a <tt>MPIVectorSpaceBase</tt>
 * subclass appropriately.
 *
 * If the <tt>getSubVector()</tt> methods are ever called with index
 * ranges outside of those of the local processor, then the default
 * implementations in <tt>Vector</tt> of all of the methods
 * (<tt>const</tt>) <tt>getSubVector()</tt>, <tt>freeSubVector()</tt>,
 * (non-<tt>const</tt>) <tt>getSubVector()</tt> and
 * <tt>commitSubVector()</tt> are called in instead.  Alternatively, a
 * subclass could provide more specialized implemenations of these
 * methods (for more efficient gather/scatter operations) if desired
 * but this should not be needed for most use cases.
 *
 * It is interesting to note that in the above use case, the explicit
 * subvector access methods call on its default implementation defined
 * in <tt>Vector</tt> (which calls on <tt>applyOp()</tt>), the
 * operator will be properly applied since the version of
 * <tt>applyOp()</tt> implemented in this class will only request
 * local vector data and hence there will only be two levels of
 * recussion for any call to an explicit subvector access method.
 * This is a truly elegant result.
 *
 * As described in the documentation for <tt>MPIVectorSpaceBase</tt>,
 * it is possible that at runtime it may be discovered that the
 * mapping of vector data to processors does not fall under this
 * design in which case the method <tt>applyOp()</tt> should be
 * overridden to handle this which will of course remove the
 * possibility of interoperability with other MPI-based vector
 * objects.
 *
 * Note that vector subclass derived from this node interface class
 * must only be directly used in SPMD mode for this to work properly.
 */
template<class Scalar>
class MPIVectorBase : virtual public Vector<Scalar> {
public:

	///
	using MultiVector<Scalar>::applyOp;
	///
	using Vector<Scalar>::applyOp;

	///
	MPIVectorBase();

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI-based vector space object for <tt>*this</tt> vector.
	 */
	virtual Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

	///
	/** Returns a non-<tt>const</tt> pointer to the beginning of the
	 * local vector data (and its stride).
	 *
	 * @param  values  [out] On output <tt>*values</tt> will point to an array of the local values.
	 * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*values)[]</tt>
	 */
	virtual void getLocalData( Scalar** values, ptrdiff_t* stride ) = 0;

	//@}

	/** @name Virtual methods with default implementations. */
	//@{

	///
	/** Returns a <tt>const</tt> pointer to the beginning of the local
	 * vector data.
	 *
	 * The default implementation performs a <tt>const_cast</tt> of <tt>this</tt> and
	 * then calls the non-<tt>const</tt> version.
	 */
	virtual void getLocalData( const Scalar** values, ptrdiff_t* stride ) const;

	//@}

	/** @name Overridden from MultiVector */
	//@{

	/// Calls <tt>mpiSpace()</tt>
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > space() const;
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
		,const int                   num_vecs
		,const Vector<Scalar>*          vecs[]
		,const int                   num_targ_vecs
		,Vector<Scalar>*                targ_vecs[]
		,RTOpPack::ReductTarget         *reduct_obj
		,const Index                    first_ele
		,const Index                    sub_dim
		,const Index                    global_offset
		) const;
	///
	void getSubVector( const Range1D& rng, RTOpPack::SubVectorT<Scalar>* sub_vec ) const;
	///
	void freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const;
	///
	void getSubVector( const Range1D& rng, RTOpPack::MutableSubVectorT<Scalar>* sub_vec );
	///
	void commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec );

	//@}

protected:

	///
	/** Subclasses should call whenever the structure of the VectorSpace changes.
	 */
	virtual void updateMpiSpace();

private:

	// ///////////////////////////////////////
	// Private data members
	
	mutable bool in_applyOp_;

	// Cached (only on vector space!)
	mutable Index  globalDim_;
	mutable Index  localOffset_;
	mutable Index  localSubDim_;

	// /////////////////////////////////////
	// Private member functions

	Range1D validateRange( const Range1D& rng_in ) const;

}; // end class MPIVectorBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_BASE_DECL_HPP
