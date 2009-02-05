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

#ifndef THYRA_SPMD_VECTOR_SPACE_BASE_HPP
#define THYRA_SPMD_VECTOR_SPACE_BASE_HPP

#include "Thyra_SpmdVectorSpaceDefaultBaseDecl.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_DefaultSpmdVectorSpaceFactory.hpp"
#include "Thyra_SpmdVectorSpaceUtilities.hpp"

namespace Thyra {


template<class Scalar>
SpmdVectorSpaceDefaultBase<Scalar>::SpmdVectorSpaceDefaultBase()
  :mapCode_(-1), defaultLocalOffset_(-1), defaultGlobalDim_(-1), localSubDim_(-1)
{}


// Virtual methods with default implementations


template<class Scalar>
Index SpmdVectorSpaceDefaultBase<Scalar>::localOffset() const
{
  return defaultLocalOffset_;
}


template<class Scalar>
Index SpmdVectorSpaceDefaultBase<Scalar>::mapCode() const
{
  return mapCode_;
}


template<class Scalar>
std::string SpmdVectorSpaceDefaultBase<Scalar>::description() const
{
  using Teuchos::RCP; using Teuchos::Comm; using Teuchos::null;
  std::ostringstream ostr;
  ostr << Teuchos::typeName(*this) << "{";
  ostr << "globalDim="<<this->dim();
  ostr << ",localSubDim="<<this->localSubDim();
  ostr << ",localOffset="<<this->localOffset();
  ostr << ",comm=";
  RCP<const Comm<Index> > comm;
  if ( (comm=this->getComm())!=null ) {
    ostr << comm->description();
  }
  else {
    ostr << "NULL";
  }
  ostr << "}";
  return ostr.str();
}


// Overridden from VectorSpaceBase

template<class Scalar>
Index SpmdVectorSpaceDefaultBase<Scalar>::dim() const
{
  return defaultGlobalDim_;
}

template<class Scalar>
Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> >
SpmdVectorSpaceDefaultBase<Scalar>::smallVecSpcFcty() const
{
  return smallVecSpcFcty_;
}

template<class Scalar>
bool SpmdVectorSpaceDefaultBase<Scalar>::isCompatible(
  const VectorSpaceBase<Scalar>& vecSpc
  ) const
{
  if( this->hasInCoreView() && vecSpc.hasInCoreView() )
    return this->dim() == vecSpc.dim();
  const SpmdVectorSpaceBase<Scalar>
    *mpiVecSpc = dynamic_cast<const SpmdVectorSpaceBase<Scalar>*>(&vecSpc);
  if(mpiVecSpc)
    return mapCode() == mpiVecSpc->mapCode();
  return false;
}

// protected

template<class Scalar>
void SpmdVectorSpaceDefaultBase<Scalar>::updateState( const Index globalDim )
{
  localSubDim_ = this->localSubDim(); 
  const Teuchos::RCP<const Teuchos::Comm<Index> >
    comm = this->getComm();
  if( localSubDim_ >= 0 ) {
    int numProc = 1;
    int procRank = 0;
    if( comm.get() ) {
      numProc = comm->getSize();
      procRank = comm->getRank();
    }
    if( numProc > 1 && (localSubDim_ < globalDim || globalDim < 0) ) {
      mapCode_ = SpmdVectorSpaceUtilities::computeMapCode(*comm,localSubDim_);
      defaultLocalOffset_
        = SpmdVectorSpaceUtilities::computeLocalOffset(*comm,localSubDim_);
      if( globalDim < 1 ) {
        defaultGlobalDim_
          = SpmdVectorSpaceUtilities::computeGlobalDim(*comm,localSubDim_);
      }
      else {
        defaultGlobalDim_ = globalDim;
        // ToDo: Perform global reduction to check that this is correct in
        // debug build
      }
    }
    else {
      // This is a serial or a locally-replicated parallel
      // vector space.
      mapCode_ = localSubDim_;
      defaultLocalOffset_ = 0;
      defaultGlobalDim_ = localSubDim_;
    }
  }
  else {
    mapCode_  = -1;     // Uninitialized!
    defaultLocalOffset_ = -1;
    defaultGlobalDim_ = -1;
  }
  smallVecSpcFcty_
    = Teuchos::rcp(new DefaultSpmdVectorSpaceFactory<Scalar>(comm));
}
  
} // end namespace Thyra

#endif // THYRA_SPMD_VECTOR_SPACE_BASE_HPP
