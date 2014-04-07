// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_MV_UQ_PCE_HPP
#define KOKKOS_MV_UQ_PCE_HPP

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "Kokkos_MV.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos Vector/MultiVector math functions
//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

// Utility function to see if a PCE is really a constant
template <typename Storage>
bool is_pce_constant(const Sacado::UQ::PCE<Storage>& x)
{
  typedef typename Storage::ordinal_type ordinal_type;
  typedef typename Storage::value_type value_type;

  // All size-1 expansions are constants
  const ordinal_type sz = x.size();
  if (sz == 1) return true;

  // Maybe use a tolerance????
  const value_type zero = 0;
  for (ordinal_type i=1; i<sz; ++i)
    if (x.fastAccessCoeff(i) != zero) return false;

  return true;
}

inline void raise_sacado_error(const char *msg)
{
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  cuda_abort(msg);
#else
  throw std::runtime_error(msg);
#endif
}

} // namespace Impl

// Rank-1 vector add with Sacado::UQ::PCE scalar type, constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::UQ::PCE<RS>*, RL, RD, RM>
V_Add( const Kokkos::View< Sacado::UQ::PCE<RS>*, RL, RD, RM >& r,
       const typename Sacado::UQ::PCE<XS>::value_type& av,
       const Kokkos::View< Sacado::UQ::PCE<XS>*, XL, XD, XM >& x,
       const typename Sacado::UQ::PCE<XS>::value_type& bv,
       const Kokkos::View< Sacado::UQ::PCE<YS>*, YL, YD, YM >& y,
       int n = -1)
{
  typedef Kokkos::View< Sacado::UQ::PCE<RS>*, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::UQ::PCE<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::UQ::PCE<YS>*, YL, YD, YM > YVector;

  typename RVector::flat_array_type r_flat = r;
  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * r.sacado_size();

  V_Add( r_flat, av, x_flat, bv, y_flat, n );

  return r;
}

// Rank-1 vector add with Sacado::UQ::PCE scalar type, non-constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::UQ::PCE<RS>*, RL, RD, RM>
V_Add( const Kokkos::View< Sacado::UQ::PCE<RS>*, RL, RD, RM >& r,
       const Sacado::UQ::PCE<XS>& av,
       const Kokkos::View< Sacado::UQ::PCE<XS>*, XL, XD, XM >& x,
       const Sacado::UQ::PCE<XS>& bv,
       const Kokkos::View< Sacado::UQ::PCE<YS>*, YL, YD, YM >& y,
       int n = -1)
{
  typedef Kokkos::View< Sacado::UQ::PCE<RS>*, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::UQ::PCE<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::UQ::PCE<YS>*, YL, YD, YM > YVector;

  if (Impl::is_pce_constant(av) && Impl::is_pce_constant(bv)) {
   return V_Add( r, av.fastAccessCoeff(0), x, bv.fastAccessCoeff(0), y, n );
  }
  else {
    Impl::raise_sacado_error("V_Add not implemented for non-constant a or b");
  }
  return r;
}

// Rank-2 vector add with Sacado::UQ::PCE scalar type, constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::UQ::PCE<RS>**, RL, RD, RM>
MV_Add( const Kokkos::View< Sacado::UQ::PCE<RS>**, RL, RD, RM >& r,
        const typename Sacado::UQ::PCE<XS>::value_type& av,
        const Kokkos::View< Sacado::UQ::PCE<XS>**, XL, XD, XM >& x,
        const typename Sacado::UQ::PCE<XS>::value_type& bv,
        const Kokkos::View< Sacado::UQ::PCE<YS>**, YL, YD, YM >& y,
        int n = -1)
{
  typedef Kokkos::View< Sacado::UQ::PCE<RS>**, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::UQ::PCE<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::UQ::PCE<YS>**, YL, YD, YM > YVector;

  typename RVector::flat_array_type r_flat = r;
  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * r.sacado_size();

  MV_Add( r_flat, av, x_flat, bv, y_flat, n );

  return r;
}

// Rank-2 vector add with Sacado::UQ::PCE scalar type, non-constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::UQ::PCE<RS>**, RL, RD, RM>
MV_Add( const Kokkos::View< Sacado::UQ::PCE<RS>**, RL, RD, RM >& r,
        const Sacado::UQ::PCE<XS>& av,
        const Kokkos::View< Sacado::UQ::PCE<XS>**, XL, XD, XM >& x,
        const Sacado::UQ::PCE<XS>& bv,
        const Kokkos::View< Sacado::UQ::PCE<YS>**, YL, YD, YM >& y,
        int n = -1)
{
  typedef Kokkos::View< Sacado::UQ::PCE<RS>**, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::UQ::PCE<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::UQ::PCE<YS>**, YL, YD, YM > YVector;

  if (Impl::is_pce_constant(av) && Impl::is_pce_constant(bv)) {
    return MV_Add( r, av.fastAccessCoeff(0), x, bv.fastAccessCoeff(0), y, n );
  }
  else {
    Impl::raise_sacado_error("MV_Add not implemented for non-constant a or b");
  }
  return r;
}

// Rank-1 dot product
template <typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
typename Details::InnerProductSpaceTraits< Sacado::UQ::PCE<XS> >::dot_type
V_Dot( const Kokkos::View< Sacado::UQ::PCE<XS>*, XL, XD, XM >& x,
       const Kokkos::View< Sacado::UQ::PCE<YS>*, YL, YD, YM >& y,
       int n = -1 )
{
  typedef Kokkos::View< Sacado::UQ::PCE<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::UQ::PCE<YS>*, YL, YD, YM > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * x.sacado_size();

  return V_Dot( x_flat, y_flat, n );
}

// Rank-2 dot product
template <typename rVector,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
MV_Dot( const rVector& r,
        const Kokkos::View< Sacado::UQ::PCE<XS>**, XL, XD, XM >& x,
        const Kokkos::View< Sacado::UQ::PCE<YS>**, YL, YD, YM >& y,
        int n = -1 )
{
  typedef Kokkos::View< Sacado::UQ::PCE<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::UQ::PCE<YS>**, YL, YD, YM > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * x.sacado_size();

  MV_Dot( r, x_flat, y_flat, n );
}

} // namespace Kokkos

#endif /* #ifndef KOKKOS_MV_UQ_PCE_HPP */
