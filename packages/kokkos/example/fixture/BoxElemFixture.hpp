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

#ifndef KOKKOS_EXAMPLE_BOXELEMFIXTURE_HPP
#define KOKKOS_EXAMPLE_BOXELEMFIXTURE_HPP

#include <utility>

#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>

#include <HexElement.hpp>
#include <BoxElemPart.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

/** \brief  Generate a distributed unstructured finite element mesh
 *          from a partitioned NX*NY*NZ box of elements.
 */
template< class Device , BoxElemPart::ElemOrder Order >
class BoxElemFixture {
public:

  enum { ElemNode = Order == BoxElemPart::ElemLinear ? 8 :
                    Order == BoxElemPart::ElemQuadratic ? 27 : 0 };

private:

  typedef Kokkos::Example::HexElement_TensorData< ElemNode > hex_data ;

  Kokkos::Example::BoxElemPart m_box_part ;

  Kokkos::View< unsigned*[3] ,        Device > m_node_grid ;
  Kokkos::View< unsigned*[ElemNode] , Device > m_elem_node ;
  Kokkos::View< unsigned*[2] ,        Device > m_recv_node ;
  Kokkos::View< unsigned*[2] ,        Device > m_send_node ;
  Kokkos::View< unsigned* ,           Device > m_send_node_id ;

  unsigned char m_elem_node_local[ ElemNode ][4] ;

public:

  typedef Kokkos::View< const unsigned * [ElemNode] , Device , Kokkos::MemoryUnmanaged > elem_node_type ;
  typedef Kokkos::View< const unsigned * [3] , Device , Kokkos::MemoryUnmanaged > node_grid_type ;
  typedef Kokkos::View< const unsigned * [2] , Device , Kokkos::MemoryUnmanaged > comm_list_type ;
  typedef Kokkos::View< const unsigned *     , Device , Kokkos::MemoryUnmanaged > send_nodeid_type ;

  KOKKOS_INLINE_FUNCTION
  unsigned node_count() const { return m_node_grid.dimension_0(); }

  KOKKOS_INLINE_FUNCTION
  unsigned elem_count() const { return m_elem_node.dimension_0(); }

  KOKKOS_INLINE_FUNCTION
  unsigned elem_node_local( unsigned inode , unsigned k ) const
    { return m_elem_node_local[inode][k] ; }

  KOKKOS_INLINE_FUNCTION
  unsigned node_grid( unsigned inode , unsigned iaxis ) const { return m_node_grid(inode,iaxis); }

  KOKKOS_INLINE_FUNCTION
  unsigned node_grid_max( unsigned iaxis ) const { return m_box_part.global_coord_max(iaxis); }

  KOKKOS_INLINE_FUNCTION
  unsigned elem_node( unsigned ielem , unsigned inode ) const { return m_elem_node(ielem,inode); }

  elem_node_type   elem_node()   const { return m_elem_node ; }
  node_grid_type   node_grid()   const { return m_node_grid ; }
  comm_list_type   recv_node()   const { return m_recv_node ; }
  comm_list_type   send_node()   const { return m_send_node ; }
  send_nodeid_type send_nodeid() const { return m_send_node_id ; }

  KOKKOS_INLINE_FUNCTION
  BoxElemFixture( const BoxElemFixture & rhs )
    : m_box_part( rhs.m_box_part )
    , m_node_grid( rhs.m_node_grid )
    , m_elem_node( rhs.m_elem_node )
    , m_recv_node( rhs.m_recv_node )
    , m_send_node( rhs.m_send_node )
    , m_send_node_id( rhs.m_send_node_id )
    {
      for ( unsigned i = 0 ; i < ElemNode ; ++i ) {
        m_elem_node_local[i][0] = rhs.m_elem_node_local[i][0] ;
        m_elem_node_local[i][1] = rhs.m_elem_node_local[i][1] ;
        m_elem_node_local[i][2] = rhs.m_elem_node_local[i][2] ;
        m_elem_node_local[i][3] = 0 ;
      }
    }

  BoxElemFixture & operator = ( const BoxElemFixture & rhs )
    {
      m_box_part      = rhs.m_box_part ;
      m_node_grid     = rhs.m_node_grid ;
      m_elem_node     = rhs.m_elem_node ;
      m_recv_node     = rhs.m_recv_node ;
      m_send_node     = rhs.m_send_node ;
      m_send_node_id  = rhs.m_send_node_id ;
     
      for ( unsigned i = 0 ; i < ElemNode ; ++i ) {
        m_elem_node_local[i][0] = rhs.m_elem_node_local[i][0] ;
        m_elem_node_local[i][1] = rhs.m_elem_node_local[i][1] ;
        m_elem_node_local[i][2] = rhs.m_elem_node_local[i][2] ;
        m_elem_node_local[i][3] = 0 ;
      }
      return *this ;
    }

  BoxElemFixture( const BoxElemPart::Decompose decompose ,
                  const unsigned global_size ,
                  const unsigned global_rank ,
                  const unsigned elem_nx ,
                  const unsigned elem_ny ,
                  const unsigned elem_nz )
  : m_box_part( Order , decompose , global_size , global_rank , elem_nx , elem_ny , elem_nz )
  , m_node_grid( "fixture_node_grid" , m_box_part.uses_node_count() )
  , m_elem_node( "fixture_elem_node" , m_box_part.uses_elem_count() )
  , m_recv_node( "fixture_recv_node" , m_box_part.recv_node_msg_count() )
  , m_send_node( "fixture_send_node" , m_box_part.send_node_msg_count() )
  , m_send_node_id( "fixture_send_node_id" , m_box_part.send_node_id_count() )
  {
    {
      const hex_data elem_data ;

      for ( unsigned i = 0 ; i < ElemNode ; ++i ) {
        m_elem_node_local[i][0] = elem_data.eval_map[i][0] ;
        m_elem_node_local[i][1] = elem_data.eval_map[i][1] ;
        m_elem_node_local[i][2] = elem_data.eval_map[i][2] ;
        m_elem_node_local[i][3] = 0 ;
      }
    }

    const size_t nwork = 
      std::max( m_recv_node.dimension_0() ,
      std::max( m_send_node.dimension_0() ,
      std::max( m_send_node_id.dimension_0() ,
      std::max( m_node_grid.dimension_0() ,
                m_elem_node.dimension_0() * m_elem_node.dimension_1() ))));

    Kokkos::parallel_for( nwork , *this );
  }


  // Initialization:

  typedef Device device_type ;

  KOKKOS_INLINE_FUNCTION
  void operator()( size_t i ) const
  {
    if ( i < m_elem_node.dimension_0() * m_elem_node.dimension_1() ) {

      const size_t ielem = i / ElemNode ;
      const size_t inode = i % ElemNode ;

      unsigned elem_coord[3] ;
      unsigned node_coord[3] ;

      m_box_part.uses_elem_coord( ielem , elem_coord );

      node_coord[0] = elem_coord[0] + m_elem_node_local[inode][0] ;
      node_coord[1] = elem_coord[1] + m_elem_node_local[inode][1] ;
      node_coord[2] = elem_coord[2] + m_elem_node_local[inode][2] ;

      m_elem_node(ielem,inode) = m_box_part.local_node_id( node_coord );
    }

    if ( i < m_node_grid.dimension_0() ) {
      unsigned node_coord[3] ;
      m_box_part.local_node_coord( i , node_coord );
      m_node_grid(i,0) = node_coord[0] ;
      m_node_grid(i,1) = node_coord[1] ;
      m_node_grid(i,2) = node_coord[2] ;
    }

    if ( i < m_recv_node.dimension_0() ) {
      m_recv_node(i,0) = m_box_part.recv_node_rank(i);
      m_recv_node(i,1) = m_box_part.recv_node_count(i);
    }

    if ( i < m_send_node.dimension_0() ) {
      m_send_node(i,0) = m_box_part.send_node_rank(i);
      m_send_node(i,1) = m_box_part.send_node_count(i);
    }

    if ( i < m_send_node_id.dimension_0() ) {
      m_send_node_id(i) = m_box_part.send_node_id(i);
    }
  }
};

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_BOXELEMFIXTURE_HPP */

