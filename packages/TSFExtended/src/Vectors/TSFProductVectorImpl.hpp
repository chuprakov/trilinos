// @HEADER
// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
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

// //////////////////////////////////////////////////////////////
// TSFProductVectorImpl.hpp


#ifndef TSFPRODUCTVECTORIMPL_HPP
#define TSFPRODUCTVECTORIMPL_HPP


#include "TSFVectorSpaceDecl.hpp"
#include "TSFProductVectorSpaceDecl.hpp"
#include "Teuchos_RefCountPtr.hpp"

using namespace TSFExtended;
using namespace Teuchos;






//========================================================================
template <class Scalar>
ProductVector<Scalar>::ProductVector(const VectorSpace<Scalar> &space)
{
  const ProductVectorSpace<Scalar>* pvs =
    dynamic_cast<const ProductVectorSpace<Scalar>* > (space.ptr().get());
  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "Space must be a ProductVectorSpace." << endl);
  space_ = space;
  numBlocks_ = pvs->numBlocks();
  vecsE_.resize(numBlocks_);
  vecsCore_.resize(numBlocks_);
  for (int i = 0; i < numBlocks_; i++)
    {
      vecsE_[i] = space_.getBlock(i).createMember();
      vecsCore_[i] = vecsE_[i].ptr();
    }
  setCore();
  isFinal_ = false;
}



//========================================================================
template <class Scalar>
void ProductVector<Scalar>::setCore()
{
  this->uninitialize();
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > vs = space_.ptr();
  const Teuchos::RefCountPtr<const Thyra::ProductVectorSpace<Scalar> > pvs
    = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpace<Scalar> >(vs);
  initialize(pvs, &(vecsCore_[0]));
}






//========================================================================
template <class Scalar>
void ProductVector<Scalar>::setBlock(int k, const Vector<Scalar> &vec)
{
  TEST_FOR_EXCEPTION(k < 0 || k > vecsE_.length(), runtime_error,
		     "index k  = " << k << " out of range");

  
  TEST_FOR_EXCEPTION(vec.space() != space_.getBlock(k), runtime_error,
		     "for k = " << k << 
		     " vec is not a member of the underlying space" << endl << 
		     "   vec.space = " << vec.space().description() << endl <<
		     "   spcae_.getBlock(k) = " << space_.getBlock(k).description() 
		     << endl);

  vecsE_[k] = vec;
  vecsCore_[k] = vec.ptr();
  setCore();
}





//========================================================================
template <class Scalar>
void ProductVector<Scalar>::finalize()
{
  isFinal_ = true;
  for (int i = 0; i < numBlocks_; i++)
    {
      TEST_FOR_EXCEPTION(vecsE_[i].ptr() == 0, runtime_error,
			 "Attempting to finalize BlockVector where block "
			 << i << " is not set" << endl);
    }
}



//========================================================================
template <class Scalar>
void ProductVector<Scalar>::testSpace(const VectorSpace<Scalar> &space, 
				      const string &method)
{
  const ProductVectorSpace<Scalar>* pvs =
    dynamic_cast<const ProductVectorSpace<Scalar>* > (space.ptr());
  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "In " << method << "space is not ProductVectorSpace"
		     << endl);
}


//========================================================================
template <class Scalar>
Vector<Scalar>& ProductVector<Scalar>::scale(const Scalar& alpha)
{
  for (int i = 0; i < numBlocks_; i++)
    {
      vecsE_[i].scale(alpha);
    }
}





//========================================================================
template <class Scalar>
string ProductVector<Scalar>::describe(int depth) const
{
  string ret = "";
  for (int i = 0; i < depth; i++)
    {
      ret.append("   ");
    }
  ret.append("ProductVector of dimension " + toString(space_.dim()) + "\n");
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.append(vecsE_[i].describe(depth+1) + "\n");
    }
  return ret;
}




//       /** 
//        * Add a scaled vector to this vector:
//        * \code
//        * this = this + alpha*x 
//        * \endcode
//        */
//       Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x)
// {
//   testSpace(x.space(), "update");

//   for (int i = 0; i < numBlocks_; i++)
//     {
//       vecsE_[i].update(alpha, x.getBlock(i));
//     }
  
// }

//       /** 
//        * Add a scaled vector to this vector times a constant:
//        * \code
//        * this = gamma*this + alpha*x 
//        * \endcode
//        */
//       Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
//                              const Scalar& gamma)
// {
//   testSpace(x.space(), "update");
//   for (int i = 0; i < numBlocks_; i++)
//     {
//       vecsE_[i].update(alpha, x.getBlock(i), gamma);
//     }
  
// }





//       /** 
//        * Add two scaled vectors to this vector times a constant:
//        * \code
//        * this = alpha*x + beta*y + gamma*this
//        * \endcode
//        */
//       Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
//                              const Scalar& beta, const Vector<Scalar>& y, 
//                              const Scalar& gamma)
// {
//   testSpace(x.space(), "update");
//   testSpace(y.space(), "update");
//   for (int i = 0; i < numBlocks_; i++)
//     {
//       vecsE_[i].update(alpha, x.getBlock(i), beta, y.getBlock(i), gamma);
//     }
  
// }





/** 
 * Copy the values of another vector into this vector
 * \code
 * this = x
 * \endcode
 */


template <class Scalar>
Vector<Scalar>& ProductVector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
{
  testSpace(x.space(), "update");
  for (int i = 0; i < numBlocks_; i++)
    {
      vecsE_[i].acceptcopyOf(x.getBlock(i));
    }

}

/** 
 * Create a new vector that is a copy of this vector 
 */
template <class Scalar>
Vector<Scalar> ProductVector<Scalar>::copy() const 
{
  Vector<Scalar> ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecsE_[i].copy());
    }
  return ret;
}

/** 
 * Element-by-element product (Matlab dot-star operator)
 */
template <class Scalar>
Vector<Scalar> ProductVector<Scalar>::dotStar(const Vector<Scalar>& other) const 
{
  testSpace(other.space(), "dotStar");
  Vector<Scalar> ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecsE_[i].dotStar(other.getBlock(i)));
    }
  return ret;
}



/** 
 * Element-by-element division (Matlab dot-slash operator)
 */
template <class Scalar>
Vector<Scalar> ProductVector<Scalar>::dotSlash(const Vector<Scalar>& other) const 
{
  testSpace(other.space(), "dotStar");
  Vector<Scalar> ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecsE_[i].dotSlash(other.getBlock(i)));
    }
  return ret;
}

/** 
 * Return element-by-element reciprocal as a new vector
 */
template <class Scalar>
Vector<Scalar> ProductVector<Scalar>::reciprocal() const 
{
  Vector<Scalar> ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecsE_[i].reciprocal());
    }
  return ret;
}

/** 
 * Return element-by-element absolute value as a new vector
 */
template <class Scalar>
Vector<Scalar> ProductVector<Scalar>::abs() const 
{
  Vector<Scalar> ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecsE_[i].abs());
    }
  return ret;
}

/** 
 * Overwrite self with element-by-element reciprocal
 */
template <class Scalar>
Vector<Scalar>& ProductVector<Scalar>::reciprocal() 
{
  for (int i = 0; i < numBlocks_; i++)
    {
      vecsE_[i].abs();
    }
  return *this;
}

/** 
 * Overwrite self with element-by-element absolute value 
 */
template <class Scalar>
Vector<Scalar>& ProductVector<Scalar>::abs() {;}

/** 
 * Set all elements to a constant value
 */
template <class Scalar>
void ProductVector<Scalar>::setToConstant(const Scalar& alpha) {;}

      
/** 
 * Take dot product with another vector
 */
template <class Scalar>
Scalar ProductVector<Scalar>::dot(const Vector<Scalar>& other) const {;}

/**
 * Compute the 1-norm of this vector
 */
template <class Scalar>
Scalar ProductVector<Scalar>::norm1() const {;}

/**
 * Compute the 2-norm of this vector
 */
template <class Scalar>
Scalar ProductVector<Scalar>::norm2() const {;}

/**
 * Compute the infinity-norm of this vector
 */
template <class Scalar>
Scalar ProductVector<Scalar>::normInf() const {;}

/**
 * Set all elements to zero 
 */
template <class Scalar>
void ProductVector<Scalar>::zero(){;}

//@}

/** \name Element loading interface */
//@{


#endif
