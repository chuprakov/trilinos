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
// TSFCoreVectorMultiVectorDecl.hpp

#ifndef TSFCORE_VECTOR_MULTI_VECTOR_DECL_HPP
#define TSFCORE_VECTOR_MULTI_VECTOR_DECL_HPP

#include "TSFCoreVectorDecl.hpp"

namespace TSFCore {

///
/** Generic adapter subclass that takes any <tt>MultiVector</tt> that
 * has only one column and turns it into a <tt>Vector</tt>.
 *
 * The purpose of this concrete subclass is to provide an
 * implementation for <tt>Vector</tt> given that a concrete
 * implementation for a <tt>MultiVector</tt> is already provided.  A
 * linear algebra library implementation should have to do almost
 * nothing to get a <tt>Vector</tt> implementation if a
 * <tt>MultiVector</tt> is already supported.  The primary purpose
 * for the use of this subclass is in the 
 */
template<class Scalar>
class VectorMultiVector : virtual public Vector<Scalar> {
public:

	/** @name Constructors/initializers/accessors */
	//@{

	/// Construct uninitialized (see the post-conditions for <tt>uninitialize()</tt>).
	VectorMultiVector();

  /// Calls <tt>initialize()</tt>.
  VectorMultiVector(
    const Teuchos::RefCountPtr<MultiVector<Scalar> > &mv
    );

  ///
  /** Initialize given a MultiVector object.
   *
   * Preconditions:<ul>
   * <li><tt>mv.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li><tt>mv->domain().get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li><tt>mv->domain()->dim()==1</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->mv().get() == mv.get()</tt>
   * <tt><tt>this->space().get() == mv->range().get()</tt>
   */
  void initialize(
    const Teuchos::RefCountPtr<MultiVector<Scalar> > &mv
    );

  ///
  /** Set to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->mv().get() == NULL</tt>
   * <tt><tt>this->space().get() == NULL</tt>
   */
  void uninitialize(
    Teuchos::RefCountPtr<MultiVector<Scalar> > *mv = NULL
    );

  /// Return smart pointer to non-const reference to underlying <tt>MultiVector</tt> object.
  Teuchos::RefCountPtr<MultiVector<Scalar> > mv();

  /// Return smart pointer to const reference to underlying <tt>MultiVector</tt> object.
  Teuchos::RefCountPtr<const MultiVector<Scalar> > mv() const;

  //@}

  /** @name Overridden from OpBase (forwarded to this->mv()) */
  //@{
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > range() const;
  ///
	bool opSupported(ETransp M_trans) const;
  //@}

  /** @name Overridden from LinearOp (forwarded to this->mv()) */
  //@{
  ///
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;
  ///
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;
  //@}

  /** @name Overridden from MultiVector (forwarded to this->mv()) */
  //@{
  ///
	Teuchos::RefCountPtr<Vector<Scalar> > col(Index j);
  ///
	Teuchos::RefCountPtr<MultiVector<Scalar> > clone_mv() const;
  ///
	Teuchos::RefCountPtr<const MultiVector<Scalar> > subView( const Range1D& col_rng ) const;
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	///
	Teuchos::RefCountPtr<const MultiVector<Scalar> > subView( const int numCols, const int cols[] ) const;
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const int numCols, const int cols[] );
  ///
	void applyOp(
		const RTOpPack::RTOpT<Scalar>   &primary_op
		,const int                      num_multi_vecs
		,const MultiVector<Scalar>*     multi_vecs[]
		,const int                      num_targ_multi_vecs
		,MultiVector<Scalar>*           targ_multi_vecs[]
		,RTOpPack::ReductTarget*        reduct_objs[]
		,const Index                    primary_first_ele
		,const Index                    primary_sub_dim
		,const Index                    primary_global_offset
		,const Index                    secondary_first_ele
		,const Index                    secondary_sub_dim
		) const;
	///
	void applyOp(
		const RTOpPack::RTOpT<Scalar>   &primary_op
		,const RTOpPack::RTOpT<Scalar>  &secondary_op
		,const int                      num_multi_vecs
		,const MultiVector<Scalar>*     multi_vecs[]
		,const int                      num_targ_multi_vecs
		,MultiVector<Scalar>*           targ_multi_vecs[]
		,RTOpPack::ReductTarget         *reduct_obj
		,const Index                    primary_first_ele
		,const Index                    primary_sub_dim
		,const Index                    primary_global_offset
		,const Index                    secondary_first_ele
		,const Index                    secondary_sub_dim
		) const;
  ///
	void getSubMultiVector(
		const Range1D                       &rowRng
		,const Range1D                      &colRng
		,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
		) const;
	///
	void freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* sub_mv ) const;
	///
	void getSubMultiVector(
		const Range1D                                &rowRng
		,const Range1D                               &colRng
		,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
		);
	///
	void commitSubMultiVector( RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv );
  //@}

	/** @name Overridden from Vector (defined in terms of this->mv()) */
	//@{
  ///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > space() const;
	///
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
	//@}

private:
	
	Teuchos::RefCountPtr<MultiVector<Scalar> > mv_;

}; // end class VectorMultiVector

// ///////////////////////////////////////////////////
// Inline members

template <class Scalar>
inline
Teuchos::RefCountPtr<MultiVector<Scalar> >
VectorMultiVector<Scalar>::mv()
{
  return mv_;
}

template <class Scalar>
inline
Teuchos::RefCountPtr<const MultiVector<Scalar> >
VectorMultiVector<Scalar>::mv() const
{
  return mv_;
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_MULTI_VECTOR_DECL_HPP
