/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_CG_SOLVE
#define KOKKOS_EXAMPLE_CG_SOLVE

#include <cmath>
#include <limits>
#include <Kokkos_View.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_MV.hpp>
#include <impl/Kokkos_Timer.hpp>

#include <Teuchos_CommHelpers.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

inline
double all_reduce( double local , const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
  double global = 0 ;
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_SUM , 1 , & local , & global );
  return global ;
}

inline
double all_reduce_max( double local , const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
  double global = 0 ;
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_MAX , 1 , & local , & global );
  return global ;
}



template< class ImportType , class SparseMatrixType , class VectorType , class TagType = void >
struct CGSolve ;


template< class ImportType , class SparseMatrixType , class VectorType >
struct CGSolve< ImportType , SparseMatrixType , VectorType ,
  typename Kokkos::Impl::enable_if<(
    Kokkos::is_view< VectorType >::value &&
    VectorType::rank == 1
  )>::type >
{
  typedef typename VectorType::value_type scalar_type ;

  size_t iteration ;
  double iter_time ;
  double norm_res ;

  CGSolve( const ImportType       & import ,
           const SparseMatrixType & A ,
           const VectorType       & b ,
           const VectorType       & x ,
           const size_t             maximum_iteration = 200 ,
           const double             tolerance = std::numeric_limits<double>::epsilon() )
    : iteration(0)
    , iter_time(0)
    , norm_res(0)
  {
    const size_t count_owned = import.count_owned ;
    const size_t count_total = import.count_owned + import.count_receive;

    // Need input vector to matvec to be owned + received
    VectorType pAll ( "cg::p" , count_total );

    VectorType p = Kokkos::subview< VectorType >( pAll , std::pair<size_t,size_t>(0,count_owned) );
    VectorType r ( "cg::r" , count_owned );
    VectorType Ap( "cg::Ap", count_owned );

    /* r = b - A * x ; */

    /* p  = x       */  Kokkos::deep_copy( p , x );
    /* import p     */  import( pAll );
    /* Ap = A * p   */  Kokkos::MV_Multiply( Ap , A , pAll );
    /* b - Ap => r  */  Kokkos::V_Add( r , 1.0 , b , -1.0 , Ap );
    /* p  = r       */  Kokkos::deep_copy( p , r );

    double old_rdot = Kokkos::Example::all_reduce( Kokkos::V_Dot( r , r ) , import.comm );

    norm_res  = sqrt( old_rdot );
    iteration = 0 ;

    Kokkos::Impl::Timer wall_clock ;

    while ( tolerance < norm_res && iteration < maximum_iteration ) {

      /* pAp_dot = dot( p , Ap = A * p ) */

      /* import p    */  import( pAll );
      /* Ap = A * p  */  Kokkos::MV_Multiply( Ap , A , pAll );

      const double pAp_dot = Kokkos::Example::all_reduce( Kokkos::V_Dot( p , Ap ) , import.comm );
      const double alpha   = old_rdot / pAp_dot ;

      /* x +=  alpha * p ;  */ Kokkos::V_Add( x ,  alpha, p  , 1.0 , x );
      /* r += -alpha * Ap ; */ Kokkos::V_Add( r , -alpha, Ap , 1.0 , r );

      const double r_dot = Kokkos::Example::all_reduce( Kokkos::V_Dot( r , r ) , import.comm );
      const double beta  = r_dot / old_rdot ;

      /* p = r + beta * p ; */ Kokkos::V_Add( p , 1.0 , r , beta , p );

      norm_res = sqrt( old_rdot = r_dot );

      ++iteration ;
    }

    iter_time = wall_clock.seconds();
  }
};

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_CG_SOLVE */


