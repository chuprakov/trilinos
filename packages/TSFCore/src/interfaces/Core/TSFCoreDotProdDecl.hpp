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
// TSFCoreDotProdDecl.hpp

#ifndef TSFCORE_DOT_PROD_DECL_HPP
#define TSFCORE_DOT_PROD_DECL_HPP

#include "TSFCoreScalarProd.hpp"

namespace TSFCore {

///
/** Implementation of a scalar product using the dot product divided by the size of the space.
 */
template<class Scalar>
class DotProd : public ScalarProd<Scalar> {
public:
	
	/** @name Overridden from ScalarProd */
	//@{

	///
	void scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const;

	//@}

}; // end class DotProd

} // end namespace TSFCore

#endif  // TSFCORE_DOT_PROD_DECL_HPP
