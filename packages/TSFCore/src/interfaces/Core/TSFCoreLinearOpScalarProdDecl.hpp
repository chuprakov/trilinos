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
// TSFCoreLinearOpScalarProdDecl.hpp

#ifndef TSFCORE_LINEAR_OP_SCALAR_PROD_DECL_HPP
#define TSFCORE_LINEAR_OP_SCALAR_PROD_DECL_HPP

#include "TSFCoreScalarProd.hpp"

namespace TSFCore {

///
/** \brief Concrete implementation of a scalar product using a
 * symmetric positive definite linear operator..
 *
 * This subclass will work with any <tt>Vector</tt> or
 * <tt>MultiVector</tt> implementation who's vector spaces are
 * compatible with the underlying linear operator object..
 *
 * \ingroup TSFCore_basic_adapter_support_grp
 */
template<class Scalar>
class LinearOpScalarProd : public ScalarProd<Scalar> {
public:

	/** @name Constructors, initializers, accessors */
	//@{

	///
	LinearOpScalarProd();

	///
	LinearOpScalarProd( const LinearOpHandle<Scalar> &op );

	///
	void initialize( const LinearOpHandle<Scalar> &op );

	///
	const LinearOpHandle<Scalar>& op() const;

	///
	void uninitialize( LinearOpHandle<Scalar> *op = NULL );

	//@}
	
	/** @name Overridden from ScalarProd */
	//@{

	///
	void scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const;

	///
	void apply(
		const EuclideanLinearOpBase<Scalar>   &M
		,const ETransp                        M_trans
		,const Vector<Scalar>                 &x
		,Vector<Scalar>                       *y
		,const Scalar                         alpha
		,const Scalar                         beta
		) const;

	///
	void apply(
		const EuclideanLinearOpBase<Scalar>   &M
		,const ETransp                        M_trans
		,const MultiVector<Scalar>            &X
		,MultiVector<Scalar>                  *Y
		,const Scalar                         alpha
		,const Scalar                         beta
		) const;

	//@}

private:

	LinearOpHandle<Scalar>  op_;

}; // end class LinearOpScalarProd

template<class Scalar>
inline
const LinearOpHandle<Scalar>& LinearOpScalarProd<Scalar>::op() const
{
	return op_;
}

} // end namespace TSFCore

#endif  // TSFCORE_LINEAR_OP_SCALAR_PROD_DECL_HPP
