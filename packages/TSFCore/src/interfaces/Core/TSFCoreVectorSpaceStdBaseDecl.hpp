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
 * This subclass defines machinery for extracting out the definition
 * of a scalar product as an object that can be replaced.  The default
 * implementation of scalar product is the dot product.  The idea is
 * that in most cases, the definition of a scalar product may be more
 * general than a specific concrete vector implementation (i.e. a
 * single scalar product may work with all serial and all MPI-based
 * vectors if, for example, it is implemented through an
 * <tt>RTOpPack::RTOpT</tt> object).  This subclass allows an
 * application code to set a specialized scalar product without having
 * marry a particular concrete vector (vector space) implementation.
 *
 * Almost every concrete <tt>VectorSpace</tt> subclass should inherit
 * from this subclass since it makes it easy for application
 * developers to to redefine the scalar product without having to
 * create a VectorSpace subclass which can have many repercussions.
 *
 * The reason that this machinery in this base subclass is separated
 * out from the base <tt>VectorSpace</tt> interface class is that,
 * first it would clutter the base interface since this machinery is
 * an implementation artifact and, second every <tt>VectorSpace</tt>
 * subclass will not utilize this machinery.
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
	/** Construct with a different scalar product.
	 *
	 * Preconditions:<ul>
	 * <li><tt>scalarProd.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getScalarProd().get() == scalarProd.get()</tt>
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
	 * <li><tt>this->getScalarProd().get() == scalarProd.get()</tt>
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
	
	/// Calls <tt>getScalarProd()->scalarProd(x,y)</tt>
	Scalar scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const;
	
	/// Calls <tt>getScalarProd()->scalarProds(X,Y,scalar_prods)</tt>
	void scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const;
	
	//@}

private:

	Teuchos::RefCountPtr<const ScalarProd<Scalar> > scalarProd_;

}; // end class VectorSpaceStdBase

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_STD_BASE_DECL_HPP
