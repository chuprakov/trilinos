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
// TSFCoreSerialVectorBaseDecl.hpp

#ifndef TSFCORE_VECTOR_SERIAL_BASE_DECL_HPP
#define TSFCORE_VECTOR_SERIAL_BASE_DECL_HPP

#include "TSFCoreVector.hpp"

namespace TSFCore {

///
/** Node subclass of serial vectors.
 *
 * This node subclass contains the an implementation of
 * <tt>applyOp()</tt> that relies on implementations of the methods
 * <tt>getSubVector()</tt>, <tt>freeSubVector()</tt> and
 * <tt>commitSubVector()</tt>.  This class also contains default
 * implementations of <tt>getSubVector()</tt>,
 * <tt>freeSubVector()</tt> and <tt>commitSubVector()</tt>.
 *
 * <b>Notes to subclass develoeprs</b>
 *
 * All that is needed to develop a concrete subclass is to define the
 * pure virtual functions <tt>space()</tt> (for which a subclass of
 * <tt>SerialVectorSpaceBase</tt> should be usee) and
 * <tt>getData()</tt>.
 */
template<class Scalar>
class SerialVectorBase : virtual public Vector<Scalar> {
public:

	///
	using Vector<Scalar>::applyOp;

	///
	SerialVectorBase();

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns a non-<tt>const</tt> pointer to the beginning of the
	 * local vector data (and its stride).
	 *
	 * @param  values  [out] On output <tt>*values</tt> will point to an array of the values.
	 * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*values)[]</tt>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>values!=NULL</tt>
	 * <li> <tt>stride!=NULL</tt>
	 * </ul>
	 */
	virtual void getData( Scalar** values, Index* stride ) = 0;

	//@}

	/** @name Virtual methods with default implementations. */
	//@{

	///
	/** Returns a <tt>const</tt> pointer to the beginning of the local
	 * vector data.
	 *
	 * @param  values  [out] On output <tt>*values</tt> will point to an array of the values.
	 * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*values)[]</tt>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>values!=NULL</tt>
	 * <li> <tt>stride!=NULL</tt>
	 * </ul>
	 *
	 * The default implementation performs a <tt>const_cast</tt> of
	 * <tt>this</tt> and then calls the non-<tt>const</tt> version of
	 * this function.  An override of this function should only be
	 * performed to save on an extra virtual function call.
	 */
	virtual void getData( const Scalar** values, Index* stride ) const;

	//@}

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
		,const int                      num_vecs
		,const Vector<Scalar>*          vecs[]
		,const int                      num_targ_vecs
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
	///
	void setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec );

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	mutable bool   in_applyOp_;

}; // end class SerialVectorBase

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_BASE_DECL_HPP
