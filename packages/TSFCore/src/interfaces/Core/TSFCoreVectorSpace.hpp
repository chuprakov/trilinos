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

// //////////////////////////////////////////////////////////////////////
// TSFCoreVectorSpace.hpp

#ifndef TSFCORE_VECTOR_SPACE_HPP
#define TSFCORE_VECTOR_SPACE_HPP

#include "TSFCoreVectorSpaceDecl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreSerialVectorSpaceFactory.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreMultiVectorCols.hpp"

namespace TSFCore {

// Helper classes

template<class Scalar>
class CopyVectorViewBack {
public:
  CopyVectorViewBack( const Vector<Scalar> *v, const RTOpPack::MutableSubVectorT<Scalar>  &raw_v )
    :v_(v), raw_v_(raw_v)
    {}
  ~CopyVectorViewBack()
    {
      RTOpPack::SubVectorT<Scalar> sv;
      v_->getSubVector(Range1D(),&sv);
      assign_entries( &raw_v_, sv );
      v_->freeSubVector(&sv);
    }
private:
  const Vector<Scalar>                       *v_;
  const RTOpPack::MutableSubVectorT<Scalar>  raw_v_;
};

template<class Scalar>
class CopyMultiVectorViewBack {
public:
  CopyMultiVectorViewBack( const MultiVector<Scalar> *mv, const RTOpPack::MutableSubMultiVectorT<Scalar>  &raw_mv )
    :mv_(mv), raw_mv_(raw_mv)
    {}
  ~CopyMultiVectorViewBack()
    {
      RTOpPack::SubMultiVectorT<Scalar> smv;
      mv_->getSubMultiVector(Range1D(),Range1D(),&smv);
      assign_entries( &raw_mv_, smv );
      mv_->freeSubMultiVector(&smv);
    }
private:
  const MultiVector<Scalar>                       *mv_;
  const RTOpPack::MutableSubMultiVectorT<Scalar>  raw_mv_;
};

// Virtual functions with default implementations

template<class Scalar>
bool VectorSpace<Scalar>::isInCore() const
{
	return false;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceFactory<Scalar> >
VectorSpace<Scalar>::smallVecSpcFcty() const
{
	return Teuchos::rcp(new SerialVectorSpaceFactory<Scalar>());
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> > 
VectorSpace<Scalar>::createMembers(int numMembers) const
{
	return Teuchos::rcp(new MultiVectorCols<Scalar> (Teuchos::rcp(this,false),this->smallVecSpcFcty()->createVecSpc(numMembers)));
}

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
VectorSpace<Scalar>::createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_v.subDim() != this->dim() );
#endif
  // Create a vector
  Teuchos::RefCountPtr<Vector<Scalar> > v = this->createMember();
  // Copy initial values in raw_v into vector
  RTOpPack::MutableSubVectorT<Scalar> sv;
  v->getSubVector(Range1D(),&sv);
  assign_entries( &sv, raw_v );
  v->commitSubVector(&sv);
  // Setup smart pointer to vector to copy view back out just before vector is destoryed
  set_extra_data(
    Teuchos::rcp(new CopyVectorViewBack<Scalar>(&*v,raw_v))
    ,"CopyVectorViewBack"
    ,&v
    ,true
    ,Teuchos::PRE_DESTROY
    );
  return v;
}

template<class Scalar>
Teuchos::RefCountPtr<const Vector<Scalar> >
VectorSpace<Scalar>::createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_v.subDim() != this->dim() );
#endif
  // Create a vector
  Teuchos::RefCountPtr<Vector<Scalar> > v = this->createMember();
  // Copy initial values in raw_v into vector
  RTOpPack::MutableSubVectorT<Scalar> sv;
  v->getSubVector(Range1D(),&sv);
  assign_entries( &sv, raw_v );
  v->commitSubVector(&sv);
  return v;
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
VectorSpace<Scalar>::createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_mv.subDim() != this->dim() );
#endif
  // Create a multi-vector
  Teuchos::RefCountPtr< MultiVector<Scalar> > mv = this->createMembers(raw_mv.numSubCols());
  // Copy initial values in raw_mv into multi-vector
  RTOpPack::MutableSubMultiVectorT<Scalar> smv;
  mv->getSubMultiVector(Range1D(),Range1D(),&smv);
  assign_entries( &smv, raw_mv );
  mv->commitSubMultiVector(&smv);
  // Setup smart pointer to multi-vector to copy view back out just before multi-vector is destoryed
  set_extra_data(
    Teuchos::rcp(new CopyMultiVectorViewBack<Scalar>(&*mv,raw_mv))
    ,"CopyMultiVectorViewBack"
    ,&mv
    ,true
    ,Teuchos::PRE_DESTROY
    );
  return mv;
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
VectorSpace<Scalar>::createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_mv.subDim() != this->dim() );
#endif
  // Create a multi-vector
  Teuchos::RefCountPtr< MultiVector<Scalar> > mv = this->createMembers(raw_mv.numSubCols());
  // Copy values in raw_mv into multi-vector
  RTOpPack::MutableSubMultiVectorT<Scalar> smv;
  mv->getSubMultiVector(Range1D(),Range1D(),&smv);
  assign_entries( &smv, raw_mv );
  mv->commitSubMultiVector(&smv);
  return mv;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
VectorSpace<Scalar>::clone() const
{
	return Teuchos::null;
}


} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_HPP
