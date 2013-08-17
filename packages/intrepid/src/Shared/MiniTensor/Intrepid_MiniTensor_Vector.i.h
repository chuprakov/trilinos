// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_Vector_i_h)
#define Intrepid_MiniTensor_Vector_i_h

namespace Intrepid
{

//
// Default constructor
//
template<typename T, Index N>
inline
Vector<T, N>::Vector() :
TensorBase<T, Store>::TensorBase(ORDER)
{
  this->dimension_ = N;
  return;
}

//
// Constructor that initializes to NaNs
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension)
{
  this->dimension_ = N;
  return;
}

///
/// Create vector from a specified value
///
template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension, ComponentValue value) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, value)
{
  this->dimension_ = N;
  return;
}

//
// Create vector from a scalar
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension, T const & s) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, s)
{
  this->dimension_ = N;
  return;
}

//
// Create vector specifying components
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(T const & s0, T const & s1)
{
  Vector<T, N> &
  self = (*this);

  self.set_dimension(2);

  self[0] = s0;
  self[1] = s1;

  return;
}

//
// Create vector specifying components
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(T const & s0, T const & s1, T const & s2)
{
  Vector<T, N> &
  self = (*this);

  self.set_dimension(3);

  self[0] = s0;
  self[1] = s1;
  self[2] = s2;

  return;
}

//
// Create vector from array
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(Index const dimension, T const * data_ptr) :
TensorBase<T, Store>::TensorBase(dimension, ORDER, data_ptr)
{
  this->dimension_ = N;
  return;
}

//
// Copy constructor
//
template<typename T, Index N>
inline
Vector<T, N>::Vector(Vector<T, N> const & v) :
TensorBase<T, Store>::TensorBase(v)
{
  return;
}

//
// Simple destructor
//
template<typename T, Index N>
inline
Vector<T, N>::~Vector()
{
  return;
}

//
// Indexing for constant vector
//
template<typename T, Index N>
inline
T const &
Vector<T, N>::operator()(Index const i) const
{
  return (*this)[i];
}

//
// Vector indexing
//
template<typename T, Index N>
inline
T &
Vector<T, N>::operator()(Index const i)
{
  return (*this)[i];
}

//
// Vector addition
//
template<typename S, typename T, Index N>
inline
Vector<typename Promote<S, T>::type, N>
operator+(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Vector<typename Promote<S, T>::type, N>
  w;

  add(u, v, w);

  return w;
}

//
// Vector subtraction
//
template<typename S, typename T, Index N>
inline
Vector<typename Promote<S, T>::type, N>
operator-(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Vector<typename Promote<S, T>::type, N>
  w;

  subtract(u, v, w);

  return w;
}

//
// Vector minus
//
template<typename T, Index N>
inline
Vector<T, N>
operator-(Vector<T, N> const & u)
{
  Vector<T, N>
  v;

  minus(u, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T, Index N>
inline
typename Promote<S, T>::type
operator*(Vector<S, N> const & u, Vector<T, N> const & v)
{
  return dot(u, v);
}

//
// Vector equality tested by components
//
template<typename T, Index N>
inline
bool
operator==(Vector<T, N> const & u, Vector<T, N> const & v)
{
  return equal(u, v);
}

//
// Vector inequality tested by components
//
template<typename T, Index N>
inline
bool
operator!=(Vector<T, N> const & u, Vector<T, N> const & v)
{
  return not_equal(u, v);
}

//
// Scalar vector product
//
template<typename S, typename T, Index N>
inline
typename lazy_disable_if< order_1234<S>, apply_vector< Promote<S,T>, N> >::type
operator*(S const & s, Vector<T, N> const & u)
{
  Vector<typename Promote<S, T>::type, N>
  v;

  scale(u, s, v);

  return v;
}

//
// Vector scalar product
//
template<typename S, typename T, Index N>
inline
typename lazy_disable_if< order_1234<S>, apply_vector< Promote<S,T>, N> >::type
operator*(Vector<T, N> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type, N>
  v;

  scale(u, s, v);

  return v;
}

//
// Vector scalar division
//
template<typename S, typename T, Index N>
inline
Vector<typename Promote<S, T>::type, N>
operator/(Vector<T, N> const & u, S const & s)
{
  Vector<typename Promote<S, T>::type, N>
  v;

  divide(u, s, v);

  return v;
}

//
// Vector dot product
//
template<typename S, typename T, Index N>
inline
typename Promote<S, T>::type
dot(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  typename Promote<S, T>::type
  s = 0.0;

  switch (dimension) {

    default:
      for (Index i = 0; i < dimension; ++i) {
        s += u(i) * v(i);
      }
      break;

    case 3:
      s = u(0) * v(0) + u(1) * v(1) + u(2) * v(2);
      break;

    case 2:
      s = u(0) * v(0) + u(1) * v(1);
      break;

  }

  return s;
}

//
// Cross product only valid for R^3.
//
template<typename S, typename T, Index N>
inline
Vector<typename Promote<S, T>::type, N>
cross(Vector<S, N> const & u, Vector<T, N> const & v)
{
  Index const
  dimension = u.get_dimension();

  assert(v.get_dimension() == dimension);

  Vector<typename Promote<S, T>::type, N>
  w(dimension);

  switch (dimension) {

    case 3:
      w(0) = u(1) * v(2) - u(2) * v(1);
      w(1) = u(2) * v(0) - u(0) * v(2);
      w(2) = u(0) * v(1) - u(1) * v(0);
      break;

    default:
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Cross product undefined for R^" << dimension;
      std::cerr << std::endl;
      exit(1);
      break;

  }

  return w;
}

//
// R^N vector 2-norm
// \return \f$ \sqrt{u \cdot u} \f$
//
template<typename T, Index N>
inline
T
norm(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
      s = std::sqrt(dot(u, u));
      break;

    case 3:
      s = std::sqrt(u(0) * u(0) + u(1) * u(1) + u(2) * u(2));
      break;

    case 2:
      s = std::sqrt(u(0) * u(0) + u(1) * u(1));
      break;

  }

  return s;
}

//
// R^N vector 2-norm square for fast distance calculations.
// \return \f$ u \cdot u \f$
//
template<typename T, Index N>
inline
T
norm_square(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
      s = dot(u, u);
      break;

    case 3:
      s = u(0) * u(0) + u(1) * u(1) + u(2) * u(2);
      break;

    case 2:
      s = u(0) * u(0) + u(1) * u(1);
      break;

  }

  return s;
}

//
// R^N vector 1-norm
// \return \f$ \sum_i |u_i| \f$
//
template<typename T, Index N>
inline
T
norm_1(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
      for (Index i = 0; i < dimension; ++i) {
        s += std::abs(u(i));
      }
      break;

    case 3:
      s = std::abs(u(0)) + std::abs(u(1)) + std::abs(u(2));
      break;

    case 2:
      s = std::abs(u(0)) + std::abs(u(1));
      break;

  }

  return s;
}

//
// R^N vector infinity-norm
// \return \f$ \max(|u_0|,...|u_i|,...|u_N|) \f$
//
template<typename T, Index N>
inline
T
norm_infinity(Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  T
  s = 0.0;

  switch (dimension) {

    default:
      for (Index i = 0; i < dimension; ++i) {
        s = std::max(s, std::abs(u(i)));
      }
      break;

    case 3:
      s = std::max(std::max(std::abs(u(0)), std::abs(u(1))), std::abs(u(2)));
      break;

    case 2:
      s = std::max(std::abs(u(0)), std::abs(u(1)));
      break;

  }

  return s;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Vector_i_h
