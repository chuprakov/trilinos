/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
 /* @HEADER@ */

#ifndef TSFVECTORSPACEIMPL_HPP
#define TSFVECTORSPACEIMPL_HPP


#include "TSFProductVectorSpaceDecl.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFDescribable.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;


 
//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::operator==(const VectorSpace<Scalar>& other) const 
{
  return isCompatible(other);  
}


//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::operator!=(const VectorSpace<Scalar>& other) const 
{
  return !(operator==(other));
}
    



//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::isCompatible(const VectorSpace<Scalar>& vecSpc) const 
{
  TEST_FOR_EXCEPTION(vecSpc.ptr().get() == 0, runtime_error,
                     "null argument in VectorSpace<Scalar>::isCompatible()");
  return ptr().get()->isCompatible(*(vecSpc.ptr().get()));
}





//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::contains(const Vector<Scalar> &vec) const
{
  return (operator==(vec.space()));
}


//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::numBlocks() const
{
  const ProductVectorSpace<Scalar>* pvs = 
    dynamic_cast<const ProductVectorSpace<Scalar>* > (this->ptr().get());
  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "Space not a ProductVectorSpace" << endl);
  return pvs->numBlocks();
}



//========================================================================
template <class Scalar>
VectorSpace<Scalar> VectorSpace<Scalar>::getBlock(const int i) const
{
  const ProductVectorSpace<Scalar>* pvs = 
    dynamic_cast<const ProductVectorSpace<Scalar>* > (this->ptr().get());
  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "Space not a ProductVectorSpace" << endl);
  return pvs->getBlock(i);
}


//========================================================================
template <class Scalar>
void VectorSpace<Scalar>::setBlock(int i, 
				   const VectorSpace<Scalar>& space)
{
  const ProductVectorSpace<Scalar>*  pvs = 
    dynamic_cast<const ProductVectorSpace<Scalar>* >  (this->ptr().get());

  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "Can't set block of vector space that is " <<
		     "not a ProductVectorSpace.");

  ProductVectorSpace<Scalar>* pvsc = const_cast<ProductVectorSpace<Scalar>*> (pvs);
  pvsc->setBlock(i, space);
}







#endif
