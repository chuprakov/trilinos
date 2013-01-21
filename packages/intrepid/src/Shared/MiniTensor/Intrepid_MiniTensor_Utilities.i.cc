// @HEADER
// ************************************************************************
//
//                    Intrepid MiniTensor Subpackage
//                 Copyright (2013) Sandia Corporation
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

#if !defined(Intrepid_MiniTensor_Utilities_i_cc)
#define Intrepid_MiniTensor_Utilities_i_cc

#include <cmath>
#include <limits>

namespace Intrepid {

  //
  // Sign function
  //
  template <typename T>
  inline
  int
  sgn(T const & s)
  {
    return (T(0) < s) - (s < T(0));
  }

  //
  // Copysign function
  //
  template<typename T>
  inline
  T
  copysign(T const & a, T const & b)
  {
    return b >= 0 ? std::abs(a) : -std::abs(a);
  }

  //
  // NaN function. Necessary to choose the proper underlying NaN
  // for non-floating-point types.
  // Assumption: non-floating-point types have a typedef that
  // determines the underlying floating-point type.
  //
  template<typename T>
  inline
  typename Sacado::ScalarType<T>::type
  not_a_number()
  {
    return
        std::numeric_limits<typename Sacado::ScalarType<T>::type>::quiet_NaN();
  }

  //
  // Machine epsilon function. Necessary to choose the proper underlying
  // machine epsilon for non-floating-point types.
  // Assumption: non-floating-point types have a typedef that
  // determines the underlying floating-point type.
  //
  template<typename T>
  inline
  typename Sacado::ScalarType<T>::type
  machine_epsilon()
  {
    return
        std::numeric_limits<typename Sacado::ScalarType<T>::type>::epsilon();
  }

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Utilities_i_cc
