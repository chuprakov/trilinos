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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <util/TaskPool.hpp>
#include <util/Parallel.hpp>
#include <util/ParallelReduce.hpp>
#include <util/ParallelInputStream.hpp>

#include <txblas/Reduction.hpp>

#include <txblas/cr4_mxv.h>
#include <txblas/CR4Matrix.hpp>

using namespace phdmesh ;

void test_fill_cr4_band(
  ParallelMachine comm ,
  const std::vector<unsigned> & partition ,
  const unsigned iband ,
  const unsigned nband ,
  const unsigned stride ,
  const double evalue ,
  std::vector<unsigned> & prefix ,
  std::vector<txblas_cr4> & matrix );

namespace {

struct FillWork {
  double   x_mag ;
  double * x_beg ;
  unsigned x_length ;
};

int task_rand_fill( void * arg , unsigned p_size , unsigned p_rank )
{
  const double r_max = RAND_MAX ;

  FillWork & w = * reinterpret_cast<FillWork*>( arg );

  const unsigned p_next = p_rank + 1 ;
  const unsigned       n = w.x_length ;
        double * const xe = w.x_beg + ( n * p_next ) / p_size ;
        double *       x  = w.x_beg + ( n * p_rank ) / p_size ;

  const double mag = w.x_mag ;

  unsigned seed = p_rank ;

  for ( ; xe != x ; ++x ) {
    const double r = rand_r( & seed );
    *x = mag * 2.0 * ( ( r / r_max ) - 0.5 );
  }
  return 0 ;
}

void timing_test_txblas(
  ParallelMachine comm , const unsigned M , const unsigned ncycle )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned m_rem   = M % p_size ;
  const unsigned m_local = M / p_size + ( p_rank < m_rem ? 1 : 0 );

  {
    std::vector<double> tmp ;
    if ( tmp.max_size() < m_local ) {
      std::cout << "timing_test_txblas cannot allocate "
                << m_local << std::endl ;
      return ;
    }
  }

  std::vector<double> values( m_local );

  double t = wall_time();

  {
    FillWork data ;
    data.x_mag = 1 ;
    data.x_beg = & values[0] ;
    data.x_length = m_local ;
    phdmesh::taskpool::run( & task_rand_fill , & data , 0 );
  }

  const double dt_fill = wall_dtime( t );

  Summation z ;

  for ( unsigned i = 0 ; i < ncycle ; ++i ) {
    txddot1( z.xdval , m_local , & values[0] );
    all_reduce( comm , ReduceSum<1>( & z ) );
  }

  double dt_txddot1 = wall_dtime( t ) / (double) ncycle ;

  all_reduce( comm , ReduceMax<1>( & dt_txddot1 ) );

  if ( p_rank == 0 ) {
    std::cout << "TIMING rand_fill[  " << M << " ] = "
              << dt_fill << " sec" << std::endl ;
    std::cout << "TIMING txddot1[ " << M << " ] = "
              << dt_txddot1 << " sec" << std::endl ;
    std::cout.flush();
  }
}

//----------------------------------------------------------------------

void test_txblas_accuracy( ParallelMachine comm ,
                           const unsigned M ,
                           const double mag )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned m_rem   = M % p_size ;
  const unsigned m_local = M / p_size + ( p_rank < m_rem ? 1 : 0 );

  std::vector<double> values( m_local );

  {
    FillWork data ;
    data.x_mag = mag ;
    data.x_beg = & values[0] ;
    data.x_length = m_local ;
    phdmesh::taskpool::run( & task_rand_fill , & data , 0 );
  }

  const double init = 1 ;
  double a = 0 ;
  Summation z ;
  if ( p_rank == 0 ) {
    a = init ;
    z += init ;
  }

  std::vector<double>::iterator i ;


  txdsum_add_array( z.xdval , m_local , & values[0] );

  for ( i = values.begin() ; i != values.end() ; ++i ) { a += *i ; }

  all_reduce( comm , ReduceSum<1>( & z ) & ReduceSum<1>( & a ) );

  for ( i = values.begin() ; i != values.end() ; ++i ) { *i = - *i ; }

  txdsum_add_array( z.xdval , m_local , & values[0] );

  for ( i = values.begin() ; i != values.end() ; ++i ) { a += *i ; }

  all_reduce( comm , ReduceSum<1>( & z ) & ReduceSum<1>( & a ) );


  const double z_val = z.value();

  if ( p_rank == 0 ) {
    std::cout << "ACCURACY txdsum_add_array[ N = "
              << M << " , mag = " << mag << " ] " << std::endl
              << "  init_value[ " << init
              << " ] , txdsum_add[ " << z_val
              << " ] , simple_sum[ " << a << " ]" << std::endl ;
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_txblas_cr4_mxv( ParallelMachine comm ,
                         const unsigned N ,
                         const unsigned nband ,
                         const unsigned stride ,
                         const unsigned ncycle ,
                         const double Value )
{
  const unsigned p_rank = parallel_machine_rank( comm );

  std::vector<unsigned> partition ;

  simple_partition( comm , N , partition );

  const unsigned n_local = partition[ p_rank + 1 ] - partition[ p_rank ];

  const double d_one = 1 ;
  std::vector<double> x_vec( n_local , d_one );
  std::vector<double> y_vec( n_local );

  double * const x = & x_vec[0] ;
  double * const y = & y_vec[0] ;

  double total[3] = { 0 , 0 , 0 };

  double t = wall_time();

  std::vector<unsigned>   cr4_prefix ;
  std::vector<txblas_cr4> cr4_matrix ;

  test_fill_cr4_band( comm, partition, 1, nband, stride, Value,
                      cr4_prefix , cr4_matrix );
  total[0] = wall_dtime( t );

  CR4Matrix matrix( comm , partition , cr4_prefix , cr4_matrix );
  total[1] = wall_dtime( t );

  for ( unsigned i = 0 ; i < ncycle ; ++i ) {
    matrix.multiply( x , y );
  }
  total[2] = wall_dtime( t );

  all_reduce( comm , ReduceMax<3>( total ) );

  if ( p_rank == 0 ) {
    std::cout << "Test CR4Matrix timing: " ;
    std::cout << "fill = " << total[0] ;
    std::cout << ", construct = " << total[1] ;
    std::cout << ", mxv = " ;
    std::cout << ( total[2] / ncycle );
    std::cout << " seconds" << std::endl ;
  }

  return ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_rbcr_mxv( phdmesh::ParallelMachine comm , std::istream & is )
{
  unsigned nsize = 1000000 ;
  unsigned nband = 50 ;
  unsigned stride = 1 ;
  unsigned ncycle = 1 ;
  double   value = 1 ;
  if ( is.good() ) { is >> nsize ; }
  if ( is.good() ) { is >> nband ; }
  if ( is.good() ) { is >> stride ; }
  if ( is.good() ) { is >> ncycle ; }
  if ( is.good() ) { is >> value ; }

  test_txblas_cr4_mxv( comm , nsize , nband , stride , ncycle , value );
}

void test_accuracy( phdmesh::ParallelMachine comm , std::istream & is )
{
  unsigned num = 100000000 ;
  double mag = 1e10 ;
  if ( is.good() ) { is >> num ; }
  if ( is.good() ) { is >> mag ; }
  test_txblas_accuracy( comm , num , mag );
}

void test_timing( phdmesh::ParallelMachine comm , std::istream & is )
{
  unsigned num = 100000000 ;
  unsigned ncycle = 10 ;
  if ( is.good() ) { is >> num ; }
  if ( is.good() ) { is >> ncycle ; }
  timing_test_txblas( comm , num , ncycle );
}

//----------------------------------------------------------------------


