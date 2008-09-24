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
 */

#include <util/ParallelReduce.hpp>

#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace phdmesh {

#if defined( HAVE_MPI )

//----------------------------------------------------------------------

void all_write_string( ParallelMachine arg_comm ,
                       std::ostream &  arg_root_os ,
                       const std::string & arg_msg )
{
  const int i_zero = 0 ;
  const int p_root = 0 ;
  const unsigned p_size = parallel_machine_size( arg_comm );
  const unsigned p_rank = parallel_machine_rank( arg_comm );
  
  int result ;

  // Gather the send counts on root processor

  int send_count = arg_msg.size();

  std::vector<int> recv_count( p_size , i_zero );

  int * const recv_count_ptr = & recv_count[0] ;

  result = MPI_Gather( & send_count , 1 , MPI_INT ,
                       recv_count_ptr , 1 , MPI_INT ,
                       p_root , arg_comm );

  if ( MPI_SUCCESS != result ) {
    std::ostringstream msg ;
    msg << "phdmesh::all_write FAILED: MPI_Gather = " << result ;
    throw std::runtime_error( msg.str() );
  }

  // Receive counts are only non-zero on the root processor:

  std::vector<int> recv_displ( p_size + 1 , i_zero );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    recv_displ[i+1] = recv_displ[i] + recv_count[i] ;
  }

  const unsigned recv_size = (unsigned) recv_displ[ p_size ] ;
 
  std::vector<char> buffer( recv_size );

  {
    const char * const send_ptr = arg_msg.c_str();
    char * const recv_ptr = recv_size ? & buffer[0] : (char *) NULL ;
    int * const recv_displ_ptr = & recv_displ[0] ;

    result = MPI_Gatherv( (void*) send_ptr, send_count, MPI_CHAR ,
                          recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                          p_root, arg_comm );
  }

  if ( MPI_SUCCESS != result ) {
    std::ostringstream msg ;
    msg << "phdmesh::all_write FAILED: MPI_Gatherv = " << result ;
    throw std::runtime_error( msg.str() );
  }

  if ( p_root == (int) p_rank ) {
    arg_root_os << std::endl ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( recv_count[i] ) {
        char * const ptr = & buffer[ recv_displ[i] ];
        arg_root_os.write( ptr , recv_count[i] );
        arg_root_os << std::endl ;
      }
    }
    arg_root_os.flush();
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void all_reduce_internal( ParallelMachine  arg_comm ,
                          ParallelReduceOp arg_op ,
                          void           * arg_in ,
                          void           * arg_out ,
                          unsigned         arg_len )
{
  MPI_Op mpi_op = MPI_OP_NULL ;

  MPI_Op_create( arg_op , 0 , & mpi_op );

  // The SUN was buggy when combinng an
  // MPI_Allreduce with a user defined operator,
  // use reduce/broadcast instead.
/*
  const int result = 
    MPI_Allreduce(arg_in,arg_out,arg_len,MPI_BYTE,mpi_op,arg_comm);
*/

  const int result_reduce =
    MPI_Reduce(arg_in,arg_out,arg_len,MPI_BYTE,mpi_op,0,arg_comm);

  const int result_bcast =
    MPI_Bcast(arg_out,arg_len,MPI_BYTE,0,arg_comm);

  MPI_Op_free( & mpi_op );

  if ( MPI_SUCCESS != result_reduce || MPI_SUCCESS != result_bcast ) {
    std::ostringstream msg ;
    msg << "phdmesh::all_reduce FAILED: MPI_Reduce = " << result_reduce
        << " MPI_Bcast = " << result_bcast ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void all_reduce_sum( ParallelMachine comm ,
                     const double * local , double * global , unsigned count )
{
  double * tmp = const_cast<double*>( local );
  MPI_Allreduce( tmp , global , count , MPI_DOUBLE , MPI_SUM , comm );
}

void all_reduce_sum( ParallelMachine comm ,
                     const float * local , float * global , unsigned count )
{
  float * tmp = const_cast<float*>( local );
  MPI_Allreduce( tmp , global , count , MPI_FLOAT , MPI_SUM , comm );
}

void all_reduce_sum( ParallelMachine comm ,
                     const int * local , int * global , unsigned count )
{
  int * tmp = const_cast<int*>( local );
  MPI_Allreduce( tmp , global , count , MPI_INT , MPI_SUM , comm );
}

void all_reduce_bor( ParallelMachine comm ,
                     const unsigned * local ,
                     unsigned * global , unsigned count )
{
  unsigned * tmp = const_cast<unsigned*>( local );
  MPI_Allreduce( tmp , global , count , MPI_UNSIGNED , MPI_BOR , comm );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

void all_write_string( ParallelMachine ,
                       std::ostream &  arg_root_os ,
                       const std::string & arg_msg )
{
  arg_root_os << arg_msg ;
}

void all_reduce_sum( ParallelMachine ,
                     const double * local , double * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

void all_reduce_sum( ParallelMachine ,
                     const float * local , float * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

void all_reduce_sum( ParallelMachine ,
                     const int * local , int * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

void all_reduce_bor( ParallelMachine ,
                     const unsigned * local ,
                     unsigned * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

//----------------------------------------------------------------------

void all_reduce_internal( ParallelMachine ,
                          ParallelReduceOp ,
                          void   * arg_in ,
                          void   * arg_out ,
                          unsigned arg_len )
{
  unsigned char * i = reinterpret_cast<unsigned char *>( arg_in );
  unsigned char * o = reinterpret_cast<unsigned char *>( arg_out );
  for ( unsigned char * const e = i + arg_len ; e != i ; ++i , ++o ) {
    *o = *i ;
  }
}

#endif

}

