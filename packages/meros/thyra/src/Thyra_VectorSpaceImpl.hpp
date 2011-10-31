/* @HEADER@ */
/* ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_VECTORSPACEIMPL_HPP
#define THYRA_VECTORSPACEIMPL_HPP

#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_VectorSpaceDecl.hpp"
#include "Teuchos_Describable.hpp"

namespace Thyra
{
  using namespace Teuchos;
  using std::ostream;

  //========================================================================
  template <class Scalar> inline
  VectorSpace<Scalar>&
  VectorSpace<Scalar>::operator=(const RCP<VectorSpaceBase<Scalar> >& vs)
  {
    this->setRcp(vs);
    return *this;
  }
 
  //========================================================================
  template <class Scalar> inline
  bool VectorSpace<Scalar>::operator==(const VectorSpace<Scalar>& other) const 
  {
    return isCompatible(other);  
  }

  //========================================================================
  template <class Scalar> inline
  bool VectorSpace<Scalar>::operator!=(const VectorSpace<Scalar>& other) const 
  {
    return !(operator==(other));
  }

  //========================================================================
  template <class Scalar> inline
  Vector<Scalar> VectorSpace<Scalar>::createMember() const 
  {
    return Thyra::createMember(this->constPtr());
  }

  //========================================================================
  template <class Scalar> inline
  RCP<MultiVectorBase<Scalar> > VectorSpace<Scalar>::createMembers(int n) const 
  {
    return Thyra::createMembers(this->constPtr(), n);
  }

  //========================================================================
  template <class Scalar> inline
  bool VectorSpace<Scalar>::isCompatible(const VectorSpace<Scalar>& vecSpc) const 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(vecSpc.constPtr().get() == 0, std::runtime_error,
                       "null argument in VectorSpace<Scalar>::isCompatible()");
    return this->constPtr().get()->isCompatible(*(vecSpc.constPtr().get()));
  }

  //========================================================================
  template <class Scalar> inline
  bool VectorSpace<Scalar>::contains(const Vector<Scalar> &vec) const
  {
    return (operator==(vec.space()));
  }

  //========================================================================
  template <class Scalar> inline
  int VectorSpace<Scalar>::numBlocks() const
  {
    const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
      dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->constPtr().get());
    if (pvs != 0)
      {
        return pvs->numBlocks();
      }
    return 1;
  }

  //========================================================================
  template <class Scalar> inline
  VectorSpace<Scalar> VectorSpace<Scalar>::getBlock(const int i) const
  {
    const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
      dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->constPtr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(pvs == 0 && numBlocks()!=1, std::runtime_error,
                       "Space not a ProductVectorSpace" << std::endl);
    if (pvs != 0)
      {
        return pvs->getBlock(i);
      }
    return *this;
  }

  // //========================================================================
  // template <class Scalar>
  // void VectorSpace<Scalar>::setBlock(int i, 
  // 				   const VectorSpace<Scalar>& space)
  // {
  //   const Thyra::ProductVectorSpace<Scalar>*  pvs = 
  //     dynamic_cast<const Thyra::ProductVectorSpace<Scalar>* >  (this->constPtr().get());

  //   TEUCHOS_TEST_FOR_EXCEPTION(pvs == 0, std::runtime_error,
  // 		     "Can't set block of std::vector space that is " <<
  // 		     "not a ProductVectorSpace.");

  //   Thyra::ProductVectorSpace<Scalar>* pvsc = const_cast<ProductVectorSpace<Scalar>*> (pvs);
  //   pvsc->setBlock(i, space);
  // }

  //========================================================================
  template <class Scalar> inline
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(const Teuchos::Array<VectorSpace<Scalar> >& spaces)
  {
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > data(spaces.size());
    for (unsigned int i=0; i<spaces.size(); i++)
      {
        data[i] = spaces[i].constPtr();
      }
    return rcp(new Thyra::DefaultProductVectorSpace<Scalar>(data.size(), &(data[0])));
  }

  //========================================================================
  template <class Scalar> inline
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1)
  {
    return productSpace<Scalar>(Teuchos::tuple(s1));
  }

  //========================================================================
  template <class Scalar> inline
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1, 
               VectorSpace<Scalar>& s2)
  {
    return productSpace<Scalar>(Teuchos::tuple(s1, s2));
  }

  //========================================================================
  template <class Scalar> inline
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1,VectorSpace<Scalar>& s2,
               VectorSpace<Scalar>& s3)
  {
    return productSpace<Scalar>(Teuchos::tuple(s1, s2, s3));
  }

  //========================================================================
  template <class Scalar> inline 
  int lowestLocallyOwnedIndex(const VectorSpace<Scalar>& s) 
  {
    TEUCHOS_TEST_FOR_EXCEPT(!isSPMD(s));
    RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace
      = rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(s.constPtr());
    return spmdSpace->localOffset();
  }

  //========================================================================
  template <class Scalar> inline
  int numLocalElements(const VectorSpace<Scalar>& s) 
  {
    TEUCHOS_TEST_FOR_EXCEPT(!isSPMD(s));
    RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace
      = rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(s.constPtr());
    return spmdSpace->localSubDim();
  }

  //========================================================================
  template <class Scalar> inline
  bool isSPMD(const VectorSpace<Scalar>& s) 
  {
    return dynamic_cast<const SpmdVectorSpaceBase<Scalar>* >(s.constPtr().get())!=0; 
  }

  //========================================================================
  template <class Scalar> inline
  bool indexIsLocal(const VectorSpace<Scalar>& s, Ordinal i) 
  {
    Range1D range(i,i);
    return s.constPtr()->hasInCoreView(range, VIEW_TYPE_DIRECT, STRIDE_TYPE_UNIT);
  }

} // namespace Thyra

#endif
