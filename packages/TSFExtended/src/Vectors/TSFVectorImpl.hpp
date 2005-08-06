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

#ifndef TSFVECTORIMPL_HPP
#define TSFVECTORIMPL_HPP


#include "TSFVectorSpaceDecl.hpp"

using namespace TSFExtended;

//===========================================================================
template <class Scalar> 
void Vector<Scalar>::setBlock(int i, const Vector<Scalar>& v)
{
  ProductVector<Scalar>* pv = 
    dynamic_cast<ProductVector<Scalar>* >(this->ptr().get());
  TEST_FOR_EXCEPTION(pv == 0, runtime_error,
		     "vector not product vector");
  pv->setBlock(i, v);
}  



//===========================================================================
template <class Scalar> 
Vector<Scalar> Vector<Scalar>::getBlock(int i) const
{
  ProductVector<Scalar>* pv = 
    dynamic_cast <ProductVector<Scalar>* >(this->ptr().get());
  TEST_FOR_EXCEPTION(pv == 0, runtime_error,
		     "vector not product vector");
  return pv->getBlock(i);
}

  


//===========================================================================
template <class Scalar> inline 
const AccessibleVector<Scalar>* Vector<Scalar>::castToAccessible() const
{
  const AccessibleVector<Scalar>* av 
    = dynamic_cast<const AccessibleVector<Scalar>*>(this->ptr().get());
  TEST_FOR_EXCEPTION(av==0, std::runtime_error,
		     "Attempted to cast non-accessible vector "
		     << *this << " to an AccessibleVector");
  return av;
}

//===========================================================================
template <class Scalar> inline 
LoadableVector<Scalar>* Vector<Scalar>::castToLoadable()
{
  LoadableVector<Scalar>* lv 
    = dynamic_cast<LoadableVector<Scalar>*>(this->ptr().get());
  TEST_FOR_EXCEPTION(lv==0, std::runtime_error,
		     "Attempted to cast non-loadable vector "
		     << *this << " to a LoadableVector");
  return lv;
}



//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::scale(const Scalar& alpha)
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  {
    TimeMonitor t(*opTimer());
    Thyra::Vt_S(p, alpha);
  }
  return *this;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
				       const Vector<Scalar>& x)
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  const Thyra::VectorBase<Scalar>* px = x.ptr().get();
  {
    TimeMonitor t(*opTimer());
    Thyra::Vp_StV(p, alpha, *px);
  }
  return *this;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  const Thyra::VectorBase<Scalar>* px = x.ptr().get();
  {
    TimeMonitor t(*opTimer());
    if (p==0) 
      {
	Vector<Scalar> me = space().createMember();
	this->ptr() = me.ptr();
      }
    Thyra::assign(p, *px);
  }
  return *this;
}

template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::copy() const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    rtn.acceptCopyOf(*this);
  }
  return rtn;
}



//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::dotStar(const Vector<Scalar>& other) const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    Thyra::ele_wise_prod(1.0, *(this->ptr)(), *(other.ptr()), rtn.ptr().get());
  }
  return rtn;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::dotSlash(const Vector<Scalar>& other) const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    Thyra::ele_wise_divide(1.0, *(this->ptr)(), *(other.ptr()), rtn.ptr().get());
  }
  return rtn;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::abs() const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    rtn.acceptCopyOf(*this);
    rtn.abs();
  }
  return rtn;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::reciprocal() const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    rtn.acceptCopyOf(*this);
    rtn.reciprocal();
  }
  return rtn;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::abs()
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  const Thyra::VectorBase<Scalar>* px = this->ptr().get();
  {
    TimeMonitor t(*opTimer());
    Thyra::abs(p, *px);
  }
  return *this;
}
  




//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::reciprocal()
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  const Thyra::VectorBase<Scalar>* px = this->ptr().get();
  {
    TimeMonitor t(*opTimer());
    Thyra::reciprocal(p, *px);
  }
  return *this;
}

  

//===========================================================================
template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
				       const Vector<Scalar>& x, 
				       const Scalar& gamma)
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  const Thyra::VectorBase<Scalar>* px = x.ptr().get();
  {
    TimeMonitor t(*opTimer());
    Thyra::linear_combination(1, &alpha, &px, gamma, p);
  }
  return *this;
}




//===========================================================================
template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
				       const Vector<Scalar>& x, 
				       const Scalar& beta, 
				       const Vector<Scalar>& y, 
				       const Scalar& gamma)
{
  Thyra::VectorBase<Scalar>* p = this->ptr().get();
  const Thyra::VectorBase<Scalar>* px = x.ptr().get();
  const Thyra::VectorBase<Scalar>* py = y.ptr().get();
  {
    TimeMonitor t(*opTimer());
    double a[2];
    a[0] = alpha;
    a[1] = beta;
    const Thyra::VectorBase<Scalar>* vecs[2];
    vecs[0] = px;
    vecs[1] = py;
    Thyra::linear_combination(2, a, vecs, gamma, p);
  }
  return *this;
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::dot(const Vector<Scalar>& other) const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::dot(*(this->ptr)(), *(other.ptr()));
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::operator*(const Vector<Scalar>& other) const 
{
  return dot(other);
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm1() const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::norm_1(*(this->ptr)());
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2() const 
{
  TimeMonitor t(*opTimer());
  return Thyra::norm_2(*(this->ptr)());
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2(const Vector<Scalar>& weights) const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::norm_2(*(weights.ptr()) ,*(this->ptr)());
}





//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::normInf() const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::norm_inf(*(this->ptr)());
}



//===========================================================================
template <class Scalar> inline 
bool Vector<Scalar>::hasNANINF() const 
{
  double x = Thyra::sum(*(this->ptr)());
  return finite(x);
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::zero()
{
  TimeMonitor t(*opTimer());
    
  Thyra::assign(this->ptr().get(), 0.0);
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::setToConstant(const Scalar& alpha)
{
  TimeMonitor t(*opTimer());
    
  Thyra::assign(this->ptr().get(), alpha);
}
  

//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max()const
{
  TimeMonitor t(*opTimer());
  return Thyra::max(*(this->ptr)());
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max(int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar maxEl;
  Scalar* maxElP = &maxEl;
  int* indexP = &index;
  Thyra::max(*(this->ptr)(), maxElP, indexP); 
  return maxEl;
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max(const Scalar& bound, int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar maxEl;
  Scalar* maxElP = &maxEl;
  int* indexP = &index;
  Thyra::maxLessThanBound(*(this->ptr)(), bound, maxElP, indexP); 
  return maxEl;

}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min()const
{
  TimeMonitor t(*opTimer());
  return Thyra::min(*(this->ptr)());
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min(int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar minEl;
  Scalar* minElP = &minEl;
  int* indexP = &index;
  Thyra::min(*(this->ptr)(), minElP, indexP); 
  return minEl;
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min(const Scalar& bound, int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar minEl;
  Scalar* minElP = &minEl;
  int* indexP = &index;
  Thyra::minGreaterThanBound(*(this->ptr)(), bound, minElP, indexP); 
  return minEl;
}


#endif
