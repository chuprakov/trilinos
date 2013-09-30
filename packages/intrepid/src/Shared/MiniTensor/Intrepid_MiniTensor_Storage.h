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

#if !defined(Intrepid_MiniTensor_Storage_h)
#define Intrepid_MiniTensor_Storage_h

#include <boost/static_assert.hpp>

#include "Intrepid_MiniTensor_Definitions.h"

namespace Intrepid {

/// Set to constant value if not dynamic
template <Index N, Index C>
struct dimension_const {
  static Index const value = C;
};

template <Index C>
struct dimension_const<DYNAMIC, C> {
  static Index const value = DYNAMIC;
};

/// Validate dimension
template <Index D>
struct check_static {
  static Index const maximum_dimension = std::numeric_limits<Index>::digits;
  BOOST_STATIC_ASSERT_MSG(D < maximum_dimension, "Dimension too large.");
  static Index const value = D;
};

template <typename Store>
inline
void
check_dynamic(Index const dimension)
{
  assert(Store::IS_DYNAMIC == true);
  assert(dimension <= std::numeric_limits<Index>::digits);
}

/// Integer power template restricted to orders defined below
template <Index D, Index O>
struct dimension_power {
  static Index const value = 0;
};

template <Index D>
struct dimension_power<D, 1> {
  static Index const value = D;
};

template <Index D>
struct dimension_power<D, 2> {
  static Index const value = D * D;
};

template <Index D>
struct dimension_power<D, 3> {
  static Index const value = D * D * D;
};

template <Index D>
struct dimension_power<D, 4> {
  static Index const value = D * D * D * D;
};

/// Integer square for manipulations between 2nd and 4rd-order tensors.
template <Index N>
struct dimension_square {
  static Index const value = 0;
};

template <>
struct dimension_square<DYNAMIC> {
  static Index const value = DYNAMIC;
};

template <>
struct dimension_square<1> {
  static Index const value = 1;
};

template <>
struct dimension_square<2> {
  static Index const value = 4;
};

template <>
struct dimension_square<3> {
  static Index const value = 9;
};

template <>
struct dimension_square<4> {
  static Index const value = 16;
};

/// Integer square root template restricted to dimensions defined below.
/// Useful for constructing a 2nd-order tensor from a 4th-order
/// tensor with static storage.
template <Index N>
struct dimension_sqrt {
  static Index const value = 0;
};

template <>
struct dimension_sqrt<DYNAMIC> {
  static Index const value = DYNAMIC;
};

template <>
struct dimension_sqrt<1> {
  static Index const value = 1;
};

template <>
struct dimension_sqrt<4> {
  static Index const value = 2;
};

template <>
struct dimension_sqrt<9> {
  static Index const value = 3;
};

template <>
struct dimension_sqrt<16> {
  static Index const value = 4;
};

/// Manipulation of static and dynamic dimensions.
template <Index N, Index P>
struct dimension_add {
  static Index const value = N + P;
};

template <Index P>
struct dimension_add<DYNAMIC, P> {
  static Index const value = DYNAMIC;
};

template <Index N, Index P>
struct dimension_subtract {
  static Index const value = N - P;
};

template <Index P>
struct dimension_subtract<DYNAMIC, P> {
  static Index const value = DYNAMIC;
};

///
/// Base static storage class. Simple linear access memory model.
///
template<typename T, Index N>
class Storage
{
public:
  typedef T value_type;
  typedef T * pointer_type;
  typedef T & reference_type;
  typedef T const * const_pointer_type;
  typedef T const & const_reference_type;

  static
  bool const
  IS_STATIC = true;

  static
  bool const
  IS_DYNAMIC = false;

  Storage() {}

  explicit
  Storage(Index const number_entries) {resize(number_entries);}

  ~Storage() {}

  T const &
  operator[](Index const i) const
  {assert(i < N); return storage_[i];}

  T &
  operator[](Index const i)
  {assert(i < N); return storage_[i];}

  Index
  size() const {return N;}

  void
  resize(Index const number_entries) {assert(number_entries == N);}

  void
  clear() {}

  pointer_type
  get_pointer() {return &storage_[0];}

  const_pointer_type
  get_const_pointer() const {return &storage_[0];}

private:

  Storage(Storage<T, N> const & s);

  Storage<T, N> &
  operator=(Storage<T, N> const & s);

  T
  storage_[N];

};

///
/// Base dynamic storage class. Simple linear access memory model.
///
template<typename T>
class Storage<T, DYNAMIC>
{
public:
  typedef T value_type;
  typedef T * pointer_type;
  typedef T & reference_type;
  typedef T const * const_pointer_type;
  typedef T const & const_reference_type;

  static
  bool const
  IS_DYNAMIC = true;

  static
  bool const
  IS_STATIC = false;

  Storage() : storage_(NULL), size_(0) {}

  explicit
  Storage(Index const number_entries) : storage_(NULL), size_(0)
  {resize(number_entries);}

  ~Storage() {clear();}

  T const &
  operator[](Index const i) const
  {assert(i < size()); return storage_[i];}

  T &
  operator[](Index const i)
  {assert(i < size()); return storage_[i];}

  Index
  size() const
  {return size_;}

  void
  resize(Index const number_entries)
  {
    if (number_entries != size_) {
      clear(); storage_ = new T[number_entries]; size_ = number_entries;
    }
  }

  void
  clear()
  {
    if (storage_ != NULL) {
      delete [] storage_; storage_ = NULL; size_ = 0;
    }
  }

  pointer_type
  get_pointer() {return storage_;}

  const_pointer_type
  get_const_pointer() const {return storage_;}

private:

  Storage(Storage<T, DYNAMIC> const & s);

  Storage<T, DYNAMIC> &
  operator=(Storage<T, DYNAMIC> const & s);

  T *
  storage_;

  Index
  size_;
};

} // namespace Intrepid

#include "Intrepid_MiniTensor_Storage.i.h"

#endif // Intrepid_MiniTensor_Storage_h
