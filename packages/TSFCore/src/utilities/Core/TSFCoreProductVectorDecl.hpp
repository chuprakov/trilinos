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

// //////////////////////////////////////////////////////////////
// TSFCoreProductVectorDecl.hpp

#ifndef TSFCORE_PRODUCT_VECTOR_DECL_HPP
#define TSFCORE_PRODUCT_VECTOR_DECL_HPP

#include "TSFCoreProductVectorBase.hpp"

namespace TSFCore {

///
template <class Scalar> class ProductVectorSpace;

///
/** Concrete implementation of a product vector.
 *
 * Objects of this type
 *
 * The default constructor is made private to avoid accidental default
 * construction.
 */
template<class Scalar>
class ProductVector : virtual public ProductVectorBase<Scalar> {
public:

	///
	//using MultiVector<Scalar>::applyOp;

	/** @name Constructors/initializers/accessors */
	//@{

	/// Constructs to initialized (calls <tt>initialize()</tt>).
	ProductVector(
		const Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >  &productSpace
		,const Teuchos::RefCountPtr<Vector<Scalar> >                   vecs[]
		);

	///
	/** Initialize.
	 *
	 * ToDo: Finish documentation.
	 */
	void initialize(
		const Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >  &productSpace
		,const Teuchos::RefCountPtr<Vector<Scalar> >                   vecs[]
		);

	///
	/** Uninitialize.
	 *
	 * ToDo: Finish documentation.
	 */
	void uninitialize(
		Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >  *productSpace = NULL
		,Teuchos::RefCountPtr<Vector<Scalar> >                   vecs[]        = NULL
		);

	//@}

	/** @name Overridden from ProductVectorBase */
	//@{

	///
	Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > productSpace() const;
	///
	Teuchos::RefCountPtr<Vector<Scalar> > getBlock(const int k); 
	///
	Teuchos::RefCountPtr<const Vector<Scalar> > getBlock(const int k) const;

	//@}

	/** @name Overridden from Vector */
	//@{

	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > space() const;
	///
	void applyOp(
		const RTOpPack::RTOpT<Scalar>    &op
		,const int                       num_vecs
		,const Vector<Scalar>*           vecs[]
		,const int                       num_targ_vecs
		,Vector<Scalar>*                 targ_vecs[]
		,RTOpPack::ReductTarget          *reduct_obj
		,const Index                     first_ele
		,const Index                     sub_dim
		,const Index                     global_offset
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

	// //////////////////////////////
	// Private data members

	Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >       productSpace_;
	std::vector<Teuchos::RefCountPtr<Vector<Scalar> > >           vecs_;
	// cache
	int numBlocks_;


protected:

	// //////////////////////////////
	// Protected member functions
  // Added to allow TSFExtended ProductVector to derive from this.
  ProductVector();

};

} // namespace TSFCore

#endif // TSFCORE_PRODUCT_VECTOR_DECL_HPP
