/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace {

template< class > class UnitTestMDArrayView ;

template<>
class UnitTestMDArrayView< Kokkos :: KOKKOS_MACRO_DEVICE >
{
public:
  typedef Kokkos:: KOKKOS_MACRO_DEVICE device ;

  typedef Kokkos::MDArrayView< double , device > dView ;
  typedef Kokkos::MDArrayView< int ,    device > iView ;

  static std::string name()
  {
    std::string tmp ;
    tmp.append( "UnitTestMDArrayView< Kokkos::" );
    tmp.append( KOKKOS_MACRO_TO_STRING( KOKKOS_MACRO_DEVICE ) );
    tmp.append( " >" );
    return tmp ;
  }

  void error( const char * msg ) const
  {
    std::string tmp = name();
    tmp.append( msg );
    throw std::runtime_error( tmp );
  }

  UnitTestMDArrayView()
  {
    enum { dN = 1000 , iN = 2000 };

    dView dx , dy ;
    iView ix , iy ;

    if ( dx.rank() != 0 ||
         dy.rank() != 0 ||
         ix.rank() != 0 ||
         iy.rank() != 0 ) {
      error("FAILED Initialize view");
    }

    dx = Kokkos::create_labeled_mdarray< dView > ( "dx" , dN );
    ix = Kokkos::create_labeled_mdarray< iView > ( "ix" , iN );

    if ( dx.rank() != 1 || dx.dimension(0) != dN ||
         dy.rank() != 0 ||
         ix.rank() != 1 || ix.dimension(0) != iN ||
         iy.rank() != 0 ) {
      error("FAILED Allocate view");
    }

  
    dx(0) = 20 ;
    ix(0) = 10 ;
  
    dView dz = dy = dx ;
    iView iz = iy = ix ;
  
    if ( & dx(0) != & dy(0) ||
         & dx(0) != & dz(0) ||
         & ix(0) != & iy(0) ||
         & ix(0) != & iz(0) ) {
      error("FAILED Assign view");
    }

    dx = dView();
    iy = iView();
  
    if ( dx.rank() != 0 ||
         & dy(0) != & dz(0) ||
         & ix(0) != & iz(0) ||
         iy.rank() != 0 ||
         dy(0) != 20 ||
         iz(0) != 10 ) {
      error("FAILED Clear view");
    }

    dz = dy = dView();
    iz = ix = iView();

    if ( dx.rank() != 0 || dy.rank() != 0 || dz.rank() != 0 ||
         ix.rank() != 0 || iy.rank() != 0 || iz.rank() != 0 ) {
      error("FAILED Clear all view");
    }
  }
};

}

/*--------------------------------------------------------------------------*/

