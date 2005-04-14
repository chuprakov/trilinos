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
    dynamic_cast<ProductVector<Scalar>* >(ptr().get());
  TEST_FOR_EXCEPTION(pv == 0, runtime_error,
		     "vector not product vector");
  pv->setBlock(i, v);
}  



//===========================================================================
template <class Scalar> 
Vector<Scalar> Vector<Scalar>::getBlock(int i) const
{
  ProductVector<Scalar>* pv = 
    dynamic_cast <ProductVector<Scalar>* >(ptr().get());
  TEST_FOR_EXCEPTION(pv == 0, runtime_error,
		     "vector not product vector");
  return pv->getBlock(i);
}

  


//===========================================================================
template <class Scalar> inline 
const AccessibleVector<Scalar>* Vector<Scalar>::castToAccessible() const
{
  const AccessibleVector<Scalar>* av 
    = dynamic_cast<const AccessibleVector<Scalar>*>(ptr().get());
  TEST_FOR_EXCEPTION(av==0, std::runtime_error,
		     "Attempted to cast non-accessible vector "
		     << *this << " to an AccessibleVector");
  return av;
}

template <class Scalar> inline 
LoadableVector<Scalar>* Vector<Scalar>::castToLoadable()
{
  LoadableVector<Scalar>* lv 
    = dynamic_cast<LoadableVector<Scalar>*>(ptr().get());
  TEST_FOR_EXCEPTION(lv==0, std::runtime_error,
		     "Attempted to cast non-loadable vector "
		     << *this << " to a LoadableVector");
  return lv;
}

template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::scale(const Scalar& alpha)
{
  TSFCore::Vector<Scalar>* p = ptr().get();
  {
    TimeMonitor t(*opTimer());
    TSFCore::Vt_S(p, alpha);
  }
  return *this;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
				       const Vector<Scalar>& x)
{
  TSFCore::Vector<Scalar>* p = ptr().get();
  const TSFCore::Vector<Scalar>* px = x.ptr().get();
  {
    TimeMonitor t(*opTimer());
    TSFCore::Vp_StV(p, alpha, *px);
  }
  return *this;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
{
  TSFCore::Vector<Scalar>* p = ptr().get();
  const TSFCore::Vector<Scalar>* px = x.ptr().get();
  {
    TimeMonitor t(*opTimer());
    if (p==0) 
      {
	Vector<Scalar> me = space().createMember();
	ptr() = me.ptr();
      }
    TSFCore::assign(p, *px);
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
    TSFCore::ele_wise_prod(1.0, *ptr(), *(other.ptr()), rtn.ptr().get());
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
    TSFCore::ele_wise_divide(1.0, *ptr(), *(other.ptr()), rtn.ptr().get());
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
  TSFCore::Vector<Scalar>* p = ptr().get();
  const TSFCore::Vector<Scalar>* px = ptr().get();
  {
    TimeMonitor t(*opTimer());
    TSFCore::abs(p, *px);
  }
  return *this;
}
  




//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::reciprocal()
{
  TSFCore::Vector<Scalar>* p = ptr().get();
  const TSFCore::Vector<Scalar>* px = ptr().get();
  {
    TimeMonitor t(*opTimer());
    TSFCore::reciprocal(p, *px);
  }
  return *this;
}

  

//===========================================================================
template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
				       const Vector<Scalar>& x, 
				       const Scalar& gamma)
{
  TSFCore::Vector<Scalar>* p = ptr().get();
  const TSFCore::Vector<Scalar>* px = x.ptr().get();
  {
    TimeMonitor t(*opTimer());
    TSFCore::linear_combination(1, &alpha, &px, gamma, p);
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
  TSFCore::Vector<Scalar>* p = ptr().get();
  const TSFCore::Vector<Scalar>* px = x.ptr().get();
  const TSFCore::Vector<Scalar>* py = y.ptr().get();
  {
    TimeMonitor t(*opTimer());
    double a[2];
    a[0] = alpha;
    a[1] = beta;
    const TSFCore::Vector<Scalar>* vecs[2];
    vecs[0] = px;
    vecs[1] = py;
    TSFCore::linear_combination(2, a, vecs, gamma, p);
  }
  return *this;
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::dot(const Vector<Scalar>& other) const 
{
  TimeMonitor t(*opTimer());
    
  return TSFCore::dot(*ptr(), *(other.ptr()));
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
    
  return TSFCore::norm_1(*ptr());
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2() const 
{
  TimeMonitor t(*opTimer());
    
  return TSFCore::norm_2(*ptr());
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2(const Vector<Scalar>& weights) const 
{
  TimeMonitor t(*opTimer());
    
  return TSFCore::norm_2(*(weights.ptr()) ,*ptr());
}





//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::normInf() const 
{
  TimeMonitor t(*opTimer());
    
  return TSFCore::norm_inf(*ptr());
}



//===========================================================================
template <class Scalar> inline 
bool Vector<Scalar>::hasNANINF() const 
{
  double x = TSFCore::sum(*ptr());
  return finite(x);
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::zero()
{
  TimeMonitor t(*opTimer());
    
  TSFCore::assign(ptr().get(), 0.0);
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::setToConstant(const Scalar& alpha)
{
  TimeMonitor t(*opTimer());
    
  TSFCore::assign(ptr().get(), alpha);
}
  


#endif
