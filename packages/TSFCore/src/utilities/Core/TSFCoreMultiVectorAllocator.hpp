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
// TSFCoreMultiVectorAllocator.hpp

#ifndef TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP
#define TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP

#include "TSFCoreVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

///
/** Allocator class to be used with <tt>Teuchos::AbstractFactoryStd</tt> to create
 * <tt>MultiVector</tt> objects of a given size.
 */
template<class Scalar>
class MultiVectorAllocator {
public:
	///
	MultiVectorAllocator() : numMembers_(0) {}
	///
	typedef Teuchos::RefCountPtr<MultiVector<Scalar> >  ptr_t;         // required!
	///
	MultiVectorAllocator( const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vs, int numMembers )
		: vs_(vs), numMembers_(numMembers)
		{
#ifdef _DEBUG
			TEST_FOR_EXCEPTION( vs.get()==NULL, std::logic_error, "Error!" );
#endif			
		}
	///
	const ptr_t allocate() const { return vs_->createMembers(numMembers_); }  // required!
private:
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >      vs_;
	int                                                        numMembers_;
};

} // namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP
