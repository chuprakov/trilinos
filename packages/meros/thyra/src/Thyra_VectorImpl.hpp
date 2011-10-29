// @HEADER
// ***********************************************************************
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
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_IMPL_HPP
#define THYRA_VECTOR_IMPL_HPP

#include "Thyra_VectorDecl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearCombinationImpl.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_as.hpp"

namespace Thyra
{

//
// ConstVector
//

  template <class Scalar> inline
  Scalar ConstVector<Scalar>::operator[](Ordinal globalIndex) const 
  {
    ConstDetachedVectorView<Scalar> view(this->constPtr(), Range1D(0, dim(*this)-1));
    return view[globalIndex];
  }

  template <class Scalar> inline
  bool ConstVector<Scalar>::containsVector(const Thyra::VectorBase<Scalar>* vec) const
  {
    return this->constPtr().get()==vec;
  }

  template <class Scalar> inline
  void ConstVector<Scalar>::evalInto(Vector<Scalar>& acceptor) const
  {
    acceptor.acceptCopyOf(*this);
  }

  template <class Scalar> inline
  void ConstVector<Scalar>::addInto(Vector<Scalar>& acceptor, LCSign sign) const
  {
    Scalar s = convertTo<Scalar>(sign);
    Thyra::axpy(s, *this, acceptor);
  }

  template <class Scalar> inline
  std::ostream& operator<<(std::ostream& os, const ConstVector<Scalar>& v)
  {
    os << v.description() ;
    return os;
  }

//
// Vector
//

  template <class Scalar> inline
  Vector<Scalar>::Vector( const VectorSpace<Scalar> &space )
  {
    *this = createMember(space);
  }

  template <class Scalar> inline
  Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const ConstVector<Scalar>& other)
  {
    Thyra::VectorBase<Scalar>* p = this->ptr().get();
    const Thyra::VectorBase<Scalar>* px = other.constPtr().get();
    
    if (p==0) 
      {
        Vector<Scalar> me = space(other).createMember();
        this->ptr() = me.ptr();
      }
    Thyra::assign(p, *px);
    return *this;
  }

  template <class Scalar> inline
  Ordinal dim(const ConstVector<Scalar>& x) 
  {
    return x.constPtr()->space()->dim();
  }

  template <class Scalar> inline
  VectorSpace<Scalar> space(const ConstVector<Scalar>& x) 
  {
    return x.constPtr()->space();
  }

  /** \brief copy */
  THYRA_UNARY_VECTOR_OP(copy, copyInto, assign, "copy")

  //===========================================================================
  template <class Scalar> inline
  int ConstVector<Scalar>::numBlocks() const
  {
    const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
      dynamic_cast <const Thyra::ProductVectorSpaceBase<Scalar>* >(space(*this).constPtr().get());
    if (pvs==0) 
      {
        return 1;
      }
    
    return pvs->numBlocks();
  }

  //===========================================================================
  template <class Scalar> inline
  ConstVector<Scalar> ConstVector<Scalar>::getBlock(Ordinal i) const
  {
    const Thyra::ProductVectorBase<Scalar>* pv = 
      dynamic_cast <const Thyra::ProductVectorBase<Scalar>* >(this->constPtr().get());
    if (pv==0) 
      {
        TEUCHOS_TEST_FOR_EXCEPTION(i != 0, std::runtime_error,
                           "Nonzero block index " << i << " into a std::vector that is not "
                           "a product std::vector");
        return *this;
      }
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > b = pv->getVectorBlock(i);
    return b;
  }

  //===========================================================================
  template <class Scalar> inline
  Vector<Scalar> Vector<Scalar>::getBlock(int i) 
  {
    Thyra::ProductVectorBase<Scalar>* pv = 
      dynamic_cast <Thyra::ProductVectorBase<Scalar>* >(this->ptr().get());
    if (pv==0) 
      {
        TEUCHOS_TEST_FOR_EXCEPTION(i != 0, std::runtime_error,
                           "Nonzero block index " << i << " into a std::vector that is not "
                           "a product std::vector");
        return *this;
      }
    Teuchos::RCP<Thyra::VectorBase<Scalar> > b = pv->getNonconstVectorBlock(i);
    return b;
  }

  //===========================================================================
  template <class Scalar> inline
  void Vector<Scalar>::setBlock(int i, const ConstVector<Scalar>& b) 
  {
    Thyra::DefaultProductVector<Scalar>* pv = 
      dynamic_cast <Thyra::DefaultProductVector<Scalar>* >(this->ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(pv == 0, std::runtime_error,
                       "setBlock() called on a std::vector that is not a default product std::vector");
    pv->setBlock(i, b.constPtr());
    
  }

  //===========================================================================
  template <class Scalar> inline
  void Vector<Scalar>::setBlock(int i, const Vector<Scalar>& b) 
  {
    Thyra::DefaultProductVector<Scalar>* pv = 
      dynamic_cast <Thyra::DefaultProductVector<Scalar>* >(this->ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(pv == 0, std::runtime_error,
                       "setBlock() called on a std::vector that is not a default product std::vector");
    pv->setNonconstBlock(i, b.ptr());
    
  }

}

#endif // THYRA_VECTOR_IMPL_HPP
