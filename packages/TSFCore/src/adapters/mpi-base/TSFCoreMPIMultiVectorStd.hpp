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

// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIMultiVectorStd.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_STD_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_STD_HPP

// Define to make some verbose output
//#define TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT

#include "TSFCoreMPIMultiVectorStdDecl.hpp"
#include "TSFCoreVectorMultiVector.hpp"

namespace TSFCore {

// Constructors/initializers/accessors

template<class Scalar>
MPIMultiVectorStd<Scalar>::MPIMultiVectorStd()
  :leadingDim_(0)
{}

template<class Scalar>
MPIMultiVectorStd<Scalar>::MPIMultiVectorStd(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiRangeSpace
  ,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >          &domainSpace
  ,const Teuchos::RefCountPtr<Scalar>                              &values
  ,const Index                                                     leadingDim
  )
{
  initialize(mpiRangeSpace,domainSpace,values,leadingDim);
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::initialize(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiRangeSpace
  ,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >          &domainSpace
  ,const Teuchos::RefCountPtr<Scalar>                              &values
  ,const Index                                                     leadingDim
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(mpiRangeSpace.get()==NULL);
  TEST_FOR_EXCEPT(domainSpace.get()==NULL);
  TEST_FOR_EXCEPT(values.get()==NULL);
  TEST_FOR_EXCEPT(leadingDim < mpiRangeSpace->localSubDim());
#endif
  mpiRangeSpace_ = mpiRangeSpace;
  domainSpace_   = domainSpace;
  values_        = values;
  leadingDim_    = leadingDim;
  updateMpiSpace();
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    *mpiRangeSpace
  ,Teuchos::RefCountPtr<const VectorSpace<Scalar> >          *domainSpace
  ,Teuchos::RefCountPtr<Scalar>                              *values
  ,Index                                                     *leadingDim
  )
{
  if(mpiRangeSpace) *mpiRangeSpace = mpiRangeSpace_;
  if(domainSpace)   *domainSpace   = domainSpace_;
  if(values)        *values        = values_;
  if(leadingDim)    *leadingDim    = leadingDim_;

  mpiRangeSpace_  = Teuchos::null;
  domainSpace_    = Teuchos::null;
  values_         = Teuchos::null;
  leadingDim_     = 0;

  updateMpiSpace();
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
MPIMultiVectorStd<Scalar>::domain() const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::domain() const called!\n";
#endif
  return domainSpace_;
}

// Overridden from MultiVector

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
MPIMultiVectorStd<Scalar>::col(Index j)
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::col() called!\n";
#endif
  return Teuchos::rcp(new VectorMultiVector<Scalar>(subView(Range1D(j,j))));
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MPIMultiVectorStd<Scalar>::subView( const Range1D& col_rng_in )
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::subView() called!\n";
#endif
	const Range1D colRng = validateColRange(col_rng_in);
  return Teuchos::rcp(
    new MPIMultiVectorStd<Scalar>(
      mpiRangeSpace_
      ,mpiRangeSpace_->smallVecSpcFcty()->createVecSpc(colRng.size())
      ,Teuchos::rcp( (&*values_) + (colRng.lbound()-1)*leadingDim_, false )
      ,leadingDim_
      )
    );
}

// Overridden from MPIMultiVectorBase

template<class Scalar>
Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >
MPIMultiVectorStd<Scalar>::mpiSpace() const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::mpiSpace() const called!\n";
#endif
  return mpiRangeSpace_;
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::getLocalData( const Scalar **values, Index *leadingDim ) const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::getLocalData() const called!\n";
#endif
  *values     = &*values_;
  *leadingDim = leadingDim_;
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::freeLocalData( const Scalar *values ) const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::freeLocalData() called!\n";
#endif
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::getLocalData( Scalar **values, Index *leadingDim )
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::getLocalData() called!\n";
#endif
  *values     = &*values_;
  *leadingDim = leadingDim_;
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::commitLocalData( Scalar *values )
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::commitLocalData() called!\n";
#endif
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_STD_HPP
