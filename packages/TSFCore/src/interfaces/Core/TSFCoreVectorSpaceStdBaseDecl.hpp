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
// TSFCoreVectorSpaceStdBaseDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_STD_BASE_DECL_HPP
#define TSFCORE_VECTOR_SPACE_STD_BASE_DECL_HPP

#include "TSFCoreVectorSpace.hpp"

namespace TSFCore {

///
/** Node base class for <tt>VectorSpace</tt> that allows the
 * definition of scalar products to be swapped in and out.
 *
 * The idea is that almost every concrete <tt>VectorSpace</tt>
 * subclass should inherit from this subclass since it makes it easy
 * to redefine the scalar product.  However, the reason that this
 * functionality is separated out from the base <tt>VectorSpace</tt>
 * interface class is that first it clutters the interface since this
 * is an implementation artifact and second since every
 * <tt>VectorSpace</tt> subclass will not utilize the feature.
 *
 */
template<class Scalar>
class VectorSpaceStdBase : virtual public VectorSpace<Scalar> {
public:

	/** @name Constructors / initializers */
	//@{

	///
	/** Construct to use dot product as the default.
	 *
	 * Postconditions:<ul>
	 * <li><tt>dynamic_cast<const DotProd<Scalar>*>(&*this->getScalarProd()) != NULL</tt>
	 * </ul>
	 */
	VectorSpaceStdBase();

	///
	/** Constructs with a different scalar product.
	 *
	 * Preconditions:<ul>
	 * <li><tt>scalarProd.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>&*this->getScalarProd() == &*scalarProd</tt>
	 * </ul>
	 */
	VectorSpaceStdBase( const Teuchos::RefCountPtr<const ScalarProd<Scalar> > &scalarProd );

	///
	/** Set a different scalar product.
	 *
	 * This function is made virtual so that subclasses can override it
	 * and take control of what happens.  However, any override should
	 * call back on this base implementation to set the actual scalar
	 * product object.
	 *
	 * Preconditions:<ul>
	 * <li><tt>scalarProd.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>&*this->getScalarProd() == &*scalarProd</tt>
	 * </ul>
	 */
	virtual void setScalarProd( const Teuchos::RefCountPtr<const ScalarProd<Scalar> > &scalarProd );

	///
	/** Return the current scalar product.
	 */
	Teuchos::RefCountPtr<const ScalarProd<Scalar> > getScalarProd() const;

	//@}

	/** @name Overridden from VectorSpace */
	//@{

	///
	Scalar scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const;
	///
	void scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const;

	//@}

private:

	Teuchos::RefCountPtr<const ScalarProd<Scalar> > scalarProd_;

}; // end class VectorSpaceStdBase

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_STD_BASE_DECL_HPP
