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

// ///////////////////////////////////////////////////////////////
// TSFCoreVectorSpaceDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_DECL_HPP
#define TSFCORE_VECTOR_SPACE_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

///
/** Abstract interface for objects that represent a space for vectors.
 *
 * This interface acts primarily as an "Abstract Factory" interface
 * for creating <tt>Vector</tt> objects using the
 * <tt>createMember()</tt> method.  A <tt>%VectorSpace</tt> can also
 * create <tt> MultiVector</tt> objects which represent a compact
 * collection of vectors.  A secondary role for <tt>%VectorSpace</tt>
 * objects is to test for compatibility of vector spaces (and the
 * vectors and linear operators using those spaces) objects using the
 * <tt>isCompatible()</tt> method.
 *
 * A <tt>%VectorSpace</tt> object can exist independent from any
 * individual <tt>Vector</tt> (or <tt>MutiVector</tt>) object; Or, a
 * <tt>%VectorSpace</tt> object can have a lifetime that is dependent
 * on a single <tt>Vector</tt> ( or <tt>MultiVector</tt>) object.  The
 * same interface can serve both roles.
 *
 * Note that a <tt>%VectorSpace</tt> object must also be able to
 * create <tt>MultiVector</tt> objects with any number of column
 * vectors, and <tt>LinearOp::domain()</tt> gives a vector space of
 * that dimension.  An interesting side effect of this design is that
 * the creation of a multi-vector provides a way for clients to create
 * vector spaces of any arbitrary (although small usually) dimension.
 * In order to give the client the same ability without having to
 * create a full multi-vector object first, the method
 * <tt>smallVecSpcFcty()</tt> is included.  The method
 * <tt>smallVecSpcFcty()</tt> returns a <tt>VectorSpaceFactory</tt>
 * object that can create (serial) vector spaces of any small
 * dimension.
 *
 * A vector space is also where the scalar product for the space is
 * defined which is computed by the <tt>scalarProd()</tt> method.  A
 * scalar product allows the vector space to introduce scaling into
 * many different types of numerical algorithms.
 *
 * If the underlying object is not initialized, then <tt>dim()==0</tt>
 * will be true and none of the other methods should be called
 * or exceptions will be thrown.
 *
 * <b>Notes to subclass developers</b>
 *
 * A subclass is only required to override three methods:
 * <tt>dim()</tt>, <tt>isCompatible()</tt> and
 * <tt>createMember()</tt>.  Note that implementing the
 * <tt>createMember()</tt> method also entails defining a concrete
 * <tt>Vector</tt> subclass.
 *
 * If a subclass can support specialized multi-vectors, then the
 * <tt>createMembers()</tt> should be need to be overridden.  Note
 * that implementing the <tt>createMembers()</tt> also entails
 * defining a concrete <tt>MultiVector</tt> subclass.  For some types
 * of concrete <tt>MultiVector</tt> subclass implementations
 * (e.g. serial multi-vectors), the same default
 * <tt>VectorSpaceFactory</tt> typed object returned from the default
 * <tt>smallVecSpcFcty()</tt> method can be used.  However, more
 * specialized <tt>MultiVector</tt> subclasses (e.g. distributed
 * parallel) will require an override of the <tt>smallVecSpcFcty()</tt>
 * method to return a specialized type of <tt>VectorSpaceFactory</tt>
 * object.
 */
template<class Scalar>
class VectorSpace {
public:

	///
	virtual ~VectorSpace() {}

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Return the dimension of the vector space.
	 *
	 * If the underlying object is not initialized, then <tt>dim()==0</tt>
	 * will be true and none of the other methods should be called.
	 */
	virtual Index dim() const = 0;

	///
	/** Compare the compatibility of two vector spaces.
	 *
	 * If this function returns <tt>true</tt>, then vectors created
	 * from either of the vector spaces will be compatible and can be
	 * combined in vector operations.
	 *
	 * Invariants:<ul>
	 * <li> [<tt>this->isCompatible(vecSpc) == true</tt>] <tt>vecSpc.isCompatible(*this) == true</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->dim() != vecSpc.dim()</tt>] <tt>return == false</tt>
	 * </ul>
	 */
	virtual bool isCompatible( const VectorSpace<Scalar>& vecSpc ) const = 0;

	///
	/** Create a vector member from the vector space.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>return->dim() == this->dim()</tt>
	 * <li> <tt>return->space()->isCompatible(*this) == true</tt>
	 * </ul>
	 *
	 * @return Returns a smart reference counted pointer to a
	 * dynamically allocated vector object.  After construction the
	 * values in the vector <tt>*return</tt> are unspecified
	 * (uninitialized).  This allows for faster execution times.  Note
	 * that <tt>return->space().get() == this</tt> need not be true.
	 */
	virtual Teuchos::RefCountPtr< Vector<Scalar> > createMember() const = 0;

	///
	/** Return the scalar product of two vectors in the vector space.
	 *
	 * Preconditions:<ul>
	 * <li><tt>x.space()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>y.space()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li>The scalar product is returned.
	 * </ul>
	 */
	virtual Scalar scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const = 0;

	///
	/** Return the scalar product of each column in two multi-vectors in the vector space.
	 *
	 * @param  X            [in] Multi-vector.
	 * @param  Y            [in] Multi-vector.
	 * @param  scalar_prod  [out] Array (length <tt>X.domain()->dim()</tt>) containing the
	 *                      scalar products <tt>scalar_prod[j-1] = this->scalarProd(*X.col(j),*Y.col(j))</tt>,
	 *                      for <tt>j = 1 ... X.domain()->dim()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>X.range()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>Y.range()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>X.domain()->isCompatible(*Y.domain())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>scalar_prod[j-1] = this->scalarProd(*X.col(j),*Y.col(j))</tt>, for <tt>j = 1 ... X.domain()->dim()</tt>
	 * </ul>
	 */
	virtual void scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const = 0;

	//@}

	/** @name Virtual functions with default implementations */
	//@{

	///
	/** Returns if all of the vector elements are cheaply accessible
	 * on this processor.
	 *
	 * The default implementation returns <tt>false</tt>.
	 */
	virtual bool isInCore() const;

	///
	/** Return a <tt>VectorSpaceFactory</tt> object for the creation
	 * of (serial) vector spaces with a small dimension.
	 *
	 * The default implementation returns
	 * <tt>dynamic_cast<SerialVectorSpaceFactory>(return.get())!=NULL</tt>.
	 * Note that if a subclass overrides <tt>createMembers()</tt> then
	 * it may also need to override this method as well.
	 */
	virtual Teuchos::RefCountPtr< const VectorSpaceFactory<Scalar> > smallVecSpcFcty() const;

	///
	/** Create a set of vector members (a <tt>MultiVector</tt>) from the vector space.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>num_vecs >= 1</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return->range()->isCompatible(*this) == true</tt>
	 * <li> <tt>return->domain()->dim() == numMembers</tt>
	 * </ul>
	 *
	 * @return This method returns a smart reference counted pointer
	 * to a dynamically allocated multi-vector object.  After
	 * construction, the values in <tt>*return</tt> are unspecified
	 * (uninitialized).  This allows for faster execution times.  Note
	 * that <tt>return->range().get()==this</tt> does not have to be true
	 * but will be in may cases.
	 *
	 * The default implementation returns
	 * <tt>dynamic_cast<MultiVectorCols>(return.get())!=NULL</tt>.
	 */
	virtual Teuchos::RefCountPtr< MultiVector<Scalar> > createMembers(int numMembers) const;

  ///
  /** Create a vector member that is a non-<tt>const</tt> view of raw data.
   *
   * @param  raw_v  [in] On input contains pointer (i.e. <tt>raw_v.values()</tt>)
   *                to array that the returned <tt>Vector</tt> wil be a view of.
   *                The data pointed to by <tt>raw_v.values()</tt> must remain
   *                valid until the returned <tt>Vector</tt> object is destroyed.
   *
   * Preconditions:<ul>
   * <li><tt>raw_v</tt> has been initiaized to memory (i.e.
   *     <tt>raw_v.subDim()!=0 && raw_v.values()!=NULL</tt>).
   * <li><tt>raw_v</tt> is *consistent* with the local storage
   *     of a vector's data.  This precondition is purposefully vaigue since
   *     this function can be used an variety of specialized use-cases.
   * </ul>
   *
   * Posconditions:<ul>
   * <li>See <tt>createMember()</tt>
   * </ul>
   *
   * It is stated here that the client can not expect that the values
   * pointed to by <tt>raw_v.values()</tt> to be changed until
   * the smart pointer returned goes out of scope.  This is to allow a
   * default implementation that temporarily copies data into and out
   * of a <tt>Vector</tt> object using explicit vector access.
   *
   * The default implementation of this function simply calls
   * <tt>createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements and then only when
   * the vector is destroyed is the data copied out of the vector and
   * back into the elements pointed to by
   * <tt>raw_v.values()</tt>.
   */
	virtual Teuchos::RefCountPtr<Vector<Scalar> > createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const;

  ///
  /** Create a vector member that is a <tt>const</tt> view of raw data.
   *
   * @param  raw_v  [in] On input contains pointer (i.e. <tt>raw_v.values()</tt>)
   *                to array that the returned <tt>Vector</tt> wil be a view of.
   *                The data pointed to by <tt>raw_v.values()</tt> must remain
   *                valid until the returned <tt>Vector</tt> object is destroyed.
   *
   * This function works exactly the same as the version that takes
   * a <tt>RTOpPack::MutableSubVectorT</tt> object except that this
   * version takes a <tt>RTOpPack::SubVectorT</tt> and returns a
   * smart pointer to a <tt>const</tt> <tt>Vector</tt> object.
   *
   * Preconditions:<ul>
   * <li>See the <tt>RTOpPack::MutableSubVectorT</tt> version of this function.
   * </ul>
   *
   * Posconditions:<ul>
   * <li>See <tt>createMember()</tt>
   * </ul>
   *
   * The default implementation of this function simply calls
   * <tt>createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements.
   */
	virtual Teuchos::RefCountPtr<const Vector<Scalar> > createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const;

  ///
  /** Create a multi-vector member that is a view of raw data.
   *
   * @param  raw_mv  [in] On input contains pointer (i.e. <tt>raw_mv.values()</tt>)
   *                 to array that the returned <tt>MultiVector</tt> wil be a view of.
   *
   * Preconditions:<ul>
   * <li><tt>raw_mv</tt> has been initiaized to memory (i.e.
   *     <tt>raw_mv.subDim()!=0 && raw_mv.values()!=NULL</tt>).
   * <li><tt>raw_mv</tt> is *consistent* with the local storage
   *     of a vector's data.  This precondition is purposefully vaigue since
   *     this function can be used an variety of specialized use-cases.
   * </ul>
   *
   * Posconditions:<ul>
   * <li>See <tt>createMembers()</tt> where <tt>numMembers==raw_mv.numSubCols()</tt>
   * </ul>
   *
   * It is stated here that the client can not expect that the values
   * pointed to by <tt>raw_mv.values()</tt> to be changed until
   * the smart pointer returned goes out of scope.  This is to allow a
   * default implementation that temporarily copies data into and out
   * of a <tt>MultiVector</tt> object using explicit vector access.
   *
   * The default implementation of this function simply calls
   * <tt>createMembers(raw_mv.numSubCols())</tt> to create a
   * multi-vector then uses the explicit element access functions to
   * set the elements and then only when the multi-vector is destroyed
   * is the data copied out of the multi-vector and back into the
   * elements pointed to by <tt>raw_mv.values()</tt>.
   */
	virtual Teuchos::RefCountPtr<MultiVector<Scalar> > createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const;

  ///
  /** Create a multi-vector member that is a <tt>const</tt> view of raw data.
   *
   * @param  raw_mv  [in] On input contains pointer (i.e. <tt>raw_mv.values()</tt>)
   *                 to array that the returned <tt>MultiVector</tt> wil be a view of.
   *                 The data pointed to by <tt>raw_mv.values()</tt> must remain
   *                 valid until the returned <tt>MultiVector</tt> object is destroyed.
   *
   * This function works exactly the same as the version that takes
   * a <tt>RTOpPack::MutableSubMultiVectorT</tt> object except that this
   * version takes a <tt>RTOpPack::SubMultiVectorT</tt> and returns a
   * smart pointer to a <tt>const</tt> <tt>MultiVector</tt> object.
   *
   * Preconditions:<ul>
   * <li>See the <tt>RTOpPack::MutableSubMultiVectorT</tt> version of this function.
   * </ul>
   *
   * Posconditions:<ul>
   * <li>See <tt>createMember()</tt>
   * </ul>
   *
   * The default implementation of this function simply calls
   * <tt>createMembers()</tt> to create a multi-vector then uses the explicit
   * element access functions to set the elements.
   */
	virtual Teuchos::RefCountPtr<const MultiVector<Scalar> > createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const;

	///
	/** Clone this object (if supported).
	 *
	 * It is allowed for <tt>return.get()==NULL</tt> which means that this
	 * capability is optional.
	 *
	 * The default implementation returns <tt>return.get()==NULL</tt>.
	 */
	virtual Teuchos::RefCountPtr< const VectorSpace<Scalar> > clone() const;

	//@}

}; // end class VectorSpace

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_DECL_HPP
