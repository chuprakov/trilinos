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

// //////////////////////////////////////////////////////////////
// TSFProductVectorSpace.hpp





//========================================================================
template <class Scalar>
ProductVector<Scalar>::ProductVector(const Vectorspace &space)
{
  const ProductVectorSpace<Scalar>* pvs =
    dynamic_cast<const ProductVectorSpace<Scalar>* > (space.ptr());
  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "Space must be a ProductVectorSpace." << endl);
  space_ = space;
  vecs_.resize(numBlocks_);
  numBlocks_(space.numBlocks());
  isFinal_(true);
  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i] = space_.getBlock[i].createMember();
    }
}





//========================================================================
template <class Scalar>
ProductVector<Scalar>::ProductVector()
  :isFinal_(false)
{
  numBlocks_ = 0;
  vecs_.resize(0);
}



//========================================================================
template <class Scalar>
void ProductVector<Scalar>::setBlock(int k, const Vector &vec)
{
  if (vecs_.length() < k+1) vecs_.resize(k+1);
  TEST_FOR_EXCEPTION(vecs_[k].ptr() != 0, runtime_error,
		     "Block k = " << k << " already set" << endl);
  vecs_[k] = vec;
  numBlocks_ = vecs_.length();
}



//========================================================================
template <class Scalar>
Vector ProductVector<Scalar>::getBlock(int k) const
{
  TEST_FOR_EXCEPTION(k < 0 || k > numBlocks_, std:out_of_range,
		     "Value of k = " << k << " is out of range" << endl);
  TEST_FOR_EXCEPTION(vecs_[k].ptr() == 0, runtime_error,
		     "Block k = " << k << " is not set" << endl);
  return vecs_[k];
}




//========================================================================
template <class Scalar>
void ProductVector<Scalar>::finalize()
{
  isFinal_ = true;
  for (int i = 0; i < numBlocks_; i++)
    {
      TEST_FOR_EXCEPTION(vecs_[i].ptr() == 0, runtime_error,
			 "Attempting to finalize BlockVector where block "
			 << i << " is not set" << endl);
    }
}



//========================================================================
template <class Scalar>
void ProductVector<Scalar>::testSpace(const VectorSpace &space, 
				      const string &method)
{
  const ProductVectorSpace<Scalar>* pvs =
    dynamic_cast<const ProductVectorSpace<Scalar>* > (x.space().ptr());
  TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
		     "In " << method << "space is not ProductVectorSpace"
		     << endl);
}


      /** \name Math operations */
      //@{
      /** Multiply this vector by a constant scalar factor 
       * \code
       * this = alpha * this;
       * \endcode
      */
      Vector<Scalar>& scale(const Scalar& alpha)
{
  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i].scale(alpha);
    }
}

      /** 
       * Add a scaled vector to this vector:
       * \code
       * this = this + alpha*x 
       * \endcode
       */
      Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x)
{
  testSpace(x.space(), "update");

  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i].update(alpha, x.getBlock(i));
    }
  
}

      /** 
       * Add a scaled vector to this vector times a constant:
       * \code
       * this = gamma*this + alpha*x 
       * \endcode
       */
      Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
                             const Scalar& gamma)
{
  testSpace(x.space(), "update");
  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i].update(alpha, x.getBlock(i), gamma);
    }
  
}





      /** 
       * Add two scaled vectors to this vector times a constant:
       * \code
       * this = alpha*x + beta*y + gamma*this
       * \endcode
       */
      Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
                             const Scalar& beta, const Vector<Scalar>& y, 
                             const Scalar& gamma)
{
  testSpace(x.space(), "update");
  testSpace(y.space(), "update");
  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i].update(alpha, x.getBlock(i), beta, y.getBlock(i), gamma);
    }
  
}





      /** 
       * Copy the values of another vector into this vector
       * \code
       * this = x
       * \endcode
       */
      Vector<Scalar>& acceptCopyOf(const Vector<Scalar>& x)
{
  testSpace(x.space(), "update");
  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i].acceptcopyOf(x.getBlock(i));
    }

}

      /** 
       * Create a new vector that is a copy of this vector 
       */
      Vector<Scalar> copy() const 
{
  Vector<Scalar> ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecs_[i].copy());
    }
  return ret;
}

      /** 
       * Element-by-element product (Matlab dot-star operator)
       */
Vector<Scalar> dotStar(const Vector<Scalar>& other) const 
{
  testSpace(other.space(), "dotStar");
  Vector ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecs_[i].dotStar(other.getBlock(i)));
    }
  return ret;
}



      /** 
       * Element-by-element division (Matlab dot-slash operator)
       */
      Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;
{
  testSpace(other.space(), "dotStar");
  Vector ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecs_[i].dotSlash(other.getBlock(i)));
    }
  return ret;
}

      /** 
       * Return element-by-element reciprocal as a new vector
       */
Vector<Scalar> reciprocal() const 
{
  Vector ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecs_[i].reciprocal());
    }
  return ret;
}

      /** 
       * Return element-by-element absolute value as a new vector
       */
      Vector<Scalar> abs() const ;
{
  Vector ret = space_.createMember();
  for (int i = 0; i < numBlocks_; i++)
    {
      ret.setBlock(i, vecs_[i].abs());
    }
  return ret;
}

      /** 
       * Overwrite self with element-by-element reciprocal
       */
      Vector<Scalar>& reciprocal() 
{
  for (int i = 0; i < numBlocks_; i++)
    {
      vecs_[i].abs());
    }
  return *this;
}

      /** 
       * Overwrite self with element-by-element absolute value 
       */
      Vector<Scalar>& abs() ;

      /** 
       * Set all elements to a constant value
       */
      void setToConstant(const Scalar& alpha) ;

      
      /** 
       * Take dot product with another vector
       */
      Scalar dot(const Vector<Scalar>& other) const ;

      /**
       * Compute the 1-norm of this vector
       */
      Scalar norm1() const ;

      /**
       * Compute the 2-norm of this vector
       */
      Scalar norm2() const ;

      /**
       * Compute the infinity-norm of this vector
       */
      Scalar normInf() const ;

      /**
       * Set all elements to zero 
       */
      void zero();

      //@}

      /** \name Element loading interface */
      //@{
