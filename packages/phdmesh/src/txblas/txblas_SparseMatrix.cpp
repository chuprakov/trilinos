/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 * @date   August 2007
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <util/ParallelComm.hpp>

#include <txblas/SparseMatrix.hpp>
#include <txblas/sparse_mxv.h>

namespace phdmesh {

void simple_partition( ParallelMachine comm ,
                       const unsigned nglobal ,
                       std::vector<unsigned> & partition )
{
  const unsigned p_size = parallel_machine_size( comm );

  // Partition rows evenly among processors

  const unsigned min_row_per_proc = nglobal / p_size ;
  const unsigned remaining_rows   = nglobal % p_size ;

  partition.resize( p_size + 1 );

  partition[0] = 0 ;
  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    const unsigned n = min_row_per_proc + ( p < remaining_rows ? 1 : 0 );
    partition[p+1] = partition[p] + n ;
  }
}

//----------------------------------------------------------------------

namespace {

#if defined( PHDMESH_HAS_MPI )

#define PARALLEL_DATATYPE_DOUBLE   MPI_DOUBLE
#define PARALLEL_DATATYPE_UNSIGNED MPI_UNSIGNED

void all_to_all( ParallelMachine comm ,
                 ParallelDatatype type ,
                 const bool sparse ,
                 void * send_buf , const int * send_disp ,
                 void * recv_buf , const int * recv_disp )
{
  static const char method[] = "txblas_SparseMatrix.cpp(all_to_all)" ;

  std::ostringstream msg ;

  const int mpi_tag = 0 ;
  const unsigned p_rank = parallel_machine_rank( comm );
  const unsigned p_size = parallel_machine_size( comm );

  if ( sparse ) {
    int type_size = 0 ;
    MPI_Type_size( type , & type_size );

    unsigned num_recv = 0 ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      if ( p != p_rank && recv_disp[p] < recv_disp[p+1] ) { ++num_recv ; }
    }

    // Post receives for specific processors with specific sizes

    MPI_Request request_null = MPI_REQUEST_NULL ;
    std::vector<MPI_Request> request( num_recv , request_null );
    std::vector<MPI_Status>  status(  num_recv );

    int result = MPI_SUCCESS ;
    unsigned count = 0 ;

    for ( unsigned p = 0 ; result == MPI_SUCCESS && p < p_size ; ++p ) {
      const unsigned recv_size = recv_disp[p+1] - recv_disp[p] ;
      if ( p != p_rank && recv_size ) {
        void * ptr = ((char *) recv_buf) + type_size * recv_disp[p] ;
        result = MPI_Irecv( ptr , recv_size , type ,
                            p , mpi_tag , comm , & request[count] );
        ++count ;
      }
    }

    if ( MPI_SUCCESS != result ) {
      msg << method << " LOCAL[" << p_rank << "] ERROR: "
          << result << " == MPI_Irecv , " ;
    }

    //------------------------------
    // Sync to allow ready sends and for a potential error

    int local_error = MPI_SUCCESS == result ? 0 : 1 ;
    int global_error = 0 ;

    result = MPI_Allreduce( & local_error , & global_error ,
                            1 , MPI_INT , MPI_SUM , comm );

    if ( MPI_SUCCESS != result ) {
      msg << method << " GLOBAL ERROR: " << result << " == MPI_Allreduce" ;
    }
    else if ( global_error ) {
      result = MPI_ERR_UNKNOWN ;
    }
    else {
      //------------------------------
      // Ready-send the buffers, rotate the send processor
      // in a simple attempt to smooth out the communication traffic.

      for ( unsigned p = 0 ; MPI_SUCCESS == result && p < p_size ; ++p ) {
        const int p_dst = ( p + p_rank ) % p_size ;
        const int send_size = send_disp[ p_dst + 1 ] - send_disp[ p_dst ] ;
        if ( ( p_dst != (int) p_rank ) && send_size ) {
          void * ptr = ((char *) send_buf) + type_size * send_disp[ p_dst ] ;
          result = MPI_Rsend( ptr, send_size, type, p_dst, mpi_tag, comm );
        }
      }

      if ( MPI_SUCCESS != result ) {
        msg << method << " LOCAL ERROR: " << result << " == MPI_Rsend , " ;
      }
      else {
        MPI_Request * const p_request = & request[0] ;
        MPI_Status  * const p_status  = & status[0] ;

        result = MPI_Waitall( num_recv , p_request , p_status );
      }
    }
  }
  else {
    std::vector<int> send_size( p_size );
    std::vector<int> recv_size( p_size );

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      send_size[p] = send_disp[p+1] - send_disp[p] ;
      recv_size[p] = recv_disp[p+1] - recv_disp[p] ;
    }
    send_size[p_rank] = 0 ;
    recv_size[p_rank] = 0 ;

    int * sdisp = const_cast<int*>( send_disp );
    int * rdisp = const_cast<int*>( recv_disp );

    MPI_Alltoallv( send_buf , & send_size[0] , sdisp , type ,
                   recv_buf , & recv_size[0] , rdisp , type , comm );

  }
}

#else

#define PARALLEL_DATATYPE_DOUBLE   0
#define PARALLEL_DATATYPE_UNSIGNED 0

void all_to_all( ParallelMachine ,
                 ParallelDatatype ,
                 const bool ,
                 void * , const int * ,
                 void * , const int * )
{}

#endif

}

//----------------------------------------------------------------------

SparseMatrix::SparseMatrix() {}
SparseMatrix::~SparseMatrix() {}

//----------------------------------------------------------------------

void SparseMatrix::multiply(
  const double * const x ,
        double * const y ) const
{
  if ( m_comm_size <= 1 ) {
    if ( 1 < m_block_size ) {
      txblas_rbcr_mxv( m_row_size ,
                       m_block_size ,
                     & m_matrix[0] ,
                     & m_blocking[0] ,
                       x , y );
    }
    else {
      txblas_cr_mxv( m_row_size ,
                   & m_matrix[0] ,
                   & m_blocking[0] ,
                     x , y );
    }
  }
  else {

    std::vector<double> x_work( m_work_disp[ m_comm_size ] );

    {
      double * d = & x_work[0] + m_work_disp[ m_comm_rank ];
      const double * s = x ;
      const double * const e = x + m_row_size ;
      while ( e != s ) { *d++ = *s++ ; }
    }

    {
      std::vector<double> x_send( m_send_disp[ m_comm_size ] );

      for ( unsigned i = 0 ; i < m_send_map.size() ; ++i ) {
        x_send[i] = x[ m_send_map[i] ];
      }

      all_to_all( m_comm , PARALLEL_DATATYPE_DOUBLE , m_sparse ,
                  & x_send[0] , & m_send_disp[0] ,
                  & x_work[0] , & m_work_disp[0] );
    }

    if ( 1 < m_block_size ) {
      txblas_rbcr_mxv( m_row_size ,
                       m_block_size ,
                     & m_matrix[0] ,
                     & m_blocking[0] ,
                     & x_work[0] , y );
    }
    else {
      txblas_cr_mxv( m_row_size ,
                   & m_matrix[0] ,
                   & m_blocking[0] ,
                     x , y );
    }
  }
}

//----------------------------------------------------------------------

void SparseMatrix::allocate(
  ParallelMachine comm ,
  const std::vector<unsigned> & partition ,
  const unsigned local_non_zeros ,
  const unsigned local_block_size )
{
  static const char method[] = "phdmesh::SparseMatrix::allocate" ;

  m_comm = comm ;
  m_comm_size = parallel_machine_size( comm );
  m_comm_rank = parallel_machine_rank( comm );

  if ( partition.size()  != 1 + m_comm_size ) {
    std::ostringstream msg ;
    msg << method << " ERROR" ;
    msg << " comm_size = " << m_comm_size ;
    msg << " + 1  !=  partition.size() = " << partition.size() ;
    throw std::invalid_argument( msg.str() );
  }

  m_sparse     = true ;
  m_partition  = partition ;
  m_row_first  = partition[ m_comm_rank ];
  m_row_size   = partition[ m_comm_rank + 1 ] - m_row_first ;
  m_block_size = local_block_size ;

  m_matrix.resize( local_non_zeros );

  if ( 1 < m_block_size ) {
    const unsigned nblock = ( m_row_size / local_block_size ) +
                            ( m_row_size % local_block_size ? 1 : 0 );

    m_blocking.resize( nblock + 1 );
  }
  else {
    m_blocking.resize( m_row_size + 1 );
  }
}

namespace {

void ordered_insert( std::vector<unsigned> & vec , unsigned val )
{
  const std::vector<unsigned>::iterator e = vec.end();
        std::vector<unsigned>::iterator i = vec.begin();

  i = std::lower_bound( i , e , val );

  if ( e == i || val != *i ) { vec.insert( i , val ); }
}

}

void SparseMatrix::commit()
{
  static const char method[] = "phdmesh::SparseMatrix::commit" ;

  if ( 1 < m_block_size ) {
    txblas_rbcr_prep( m_row_size , m_block_size , m_matrix.size() ,
                      & m_matrix[0] , & m_blocking[0] );
  }
  else {
    txblas_cr_prep( m_row_size , m_matrix.size() ,
                    & m_matrix[0] , & m_blocking[0] );
  }

  if ( m_comm_size <= 1 ) return ;

  //------------------------------------

  const unsigned row_last = m_partition[ m_comm_rank + 1 ] - 1 ;

  m_send_disp.resize( m_comm_size + 1 );
  m_work_disp.resize( m_comm_size + 1 );

  // Generate a vector of column identifiers

  std::vector<unsigned> work_col_ident ;

  work_col_ident.reserve( m_row_size );

  {
    const std::vector<txblas_SparseMatrixEntry>::iterator j = m_matrix.end();
          std::vector<txblas_SparseMatrixEntry>::iterator b = m_matrix.begin();
          std::vector<txblas_SparseMatrixEntry>::iterator i ;

    for ( i = b ; j != i ; ++i ) { 
      if ( i->col < m_row_first ) {
        ordered_insert( work_col_ident , i->col );
      }
    }

    for ( unsigned r = 0 ; r < m_row_size ; ++r ) {
      const unsigned tmp = r + m_row_first ;
      work_col_ident.push_back( tmp );
    }

    for ( i = b ; j != i ; ++i ) { 
      if ( row_last < i->col ) {
        ordered_insert( work_col_ident , i->col );
      }
    }
  }

  //------------------------------------
  // Map column global identifiers to local work offsets

  {
    const std::vector<unsigned>::iterator b = work_col_ident.begin();
    const std::vector<unsigned>::iterator e = work_col_ident.end();

    for ( std::vector<txblas_SparseMatrixEntry>::iterator
          i = m_matrix.begin() ; i != m_matrix.end() ; ++i ) {
      const std::vector<unsigned>::iterator
        j = std::lower_bound( b, e, i->col );
      i->col = j - b ;
    }
  }

  // Displacement prefix for work vector

  {
    std::vector<unsigned>::const_iterator i = work_col_ident.begin() ;

    m_work_disp[0] = 0 ;

    for ( unsigned p = 0 ; p < m_comm_size ; ++p ) {
      const unsigned p_row_end = m_partition[p+1] ;
      unsigned count = 0 ;
      for ( ; i != work_col_ident.end() && *i < p_row_end ; ++i ) {
        ++count ;
      }

      m_work_disp[p+1] = m_work_disp[p] + count ;
    }
  }

  //------------------------------------
  // Set up communications to gather work subvector

  {
    std::vector<unsigned> send_col_size( m_comm_size );
    std::vector<unsigned> recv_col_size( m_comm_size );

    for ( unsigned p = 0 ; p < m_comm_size ; ++p ) {
      send_col_size[p] = m_work_disp[p+1] - m_work_disp[p] ;
    }
    send_col_size[ m_comm_rank ] = 0 ;

    unsigned num_msg_maximum = 0 ;

    comm_sizes( m_comm , m_comm_size / 4 , num_msg_maximum ,
                & send_col_size[0] , & recv_col_size[0] );

    m_sparse = num_msg_maximum < ( m_comm_size / 4 );

    m_send_disp[0] = 0 ;
    for ( unsigned p = 0 ; p < m_comm_size ; ++p ) {
      m_send_disp[p+1] = m_send_disp[p] + recv_col_size[p] ;
    }
  }

  const unsigned send_map_size = m_send_disp[ m_comm_size ];

  m_send_map.resize( send_map_size );

  all_to_all( m_comm , PARALLEL_DATATYPE_UNSIGNED , m_sparse ,
              & work_col_ident[0] , & m_work_disp[0],
              & m_send_map[0] ,     & m_send_disp[0] );

  for ( unsigned i = 0 ; i < send_map_size ; ++i ) {
    if ( m_send_map[i] < (int) m_row_first ||
                         (int) row_last < m_send_map[i] ) {
      std::ostringstream msg ;
      msg << method << " ERROR Received index " ;
      msg << m_send_map[i] ;
      msg << " out of range [ " ;
      msg << m_row_first ;
      msg << " : " ;
      msg << row_last ;
      msg << " ]" ;
      throw std::runtime_error( msg.str() );
    }
    m_send_map[i] -= m_row_first ;
  }
}

//----------------------------------------------------------------------

}

