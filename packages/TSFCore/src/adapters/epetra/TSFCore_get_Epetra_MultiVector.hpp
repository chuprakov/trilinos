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

// ///////////////////////////////////////////////////////
// TSFCoreget_Epetra_MultiVector.hpp

#ifndef TSFCORE_GET_EPETRA_MULTI_VECTOR_HPP
#define TSFCORE_GET_EPETRA_MULTI_VECTOR_HPP

#include "TSFCoreEpetraTypes.hpp"

namespace TSFCore {

///
/** Get a non-<tt>const</tt> <tt>Epetra_MultiVector</tt> object from a
 * non-<tt>const</tt> <tt>MultiVector</tt> object if possible.
 *
 * Preconditions:<ul>
 * <li> <tt>vs.isCompatible(*mv->range()) == true</tt>
 * </ul>
 *
 * Note: the <tt>mv</tt> object is not guaranteed to be modified until
 * the last smart pointer to the returned <tt>Epetra_MultiVector</tt>
 * object is destroyed.
 *
 * \ingroup TSFCore_adapters_Epetra_support_grp
 */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_Epetra_MultiVector(
	const EpetraVectorSpace                              &vs
	,const Teuchos::RefCountPtr<MultiVector<double> >    &mv
	);

///
/** Get a <tt>const</tt> <tt>Epetra_MultiVector</tt> object from a
 * <tt>const</tt> <tt>MultiVector</tt> object if possible.
 *
 * Preconditions:<ul>
 * <li> <tt>vs.isCompatible(*mv->range()) == true</tt>
 * </ul>
 *
 * \ingroup TSFCore_adapters_Epetra_support_grp
 */
Teuchos::RefCountPtr<const Epetra_MultiVector>
get_Epetra_MultiVector(
	const EpetraVectorSpace                                   &vs 
	,const Teuchos::RefCountPtr<const MultiVector<double> >   &mv
	);

} // namespace TSFCore

#endif // TSFCORE_GET_EPETRA_MULTI_VECTOR_HPP
