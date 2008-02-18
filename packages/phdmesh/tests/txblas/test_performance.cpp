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

#include <util/TPI.hpp>
#include <util/Parallel.hpp>
#include <util/ParallelReduce.hpp>
#include <util/ParallelInputStream.hpp>

#include <txblas/Reduction.hpp>

#include <txblas/cr4_mxv.h>
#include <txblas/CR4Matrix.hpp>

extern "C" {
void tdaxpby( TPI_ThreadPool pool ,
              unsigned n , double a , const double * x ,
                           double b , double * y );
}

using namespace phdmesh ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void test_fill_cr_band(
  ParallelMachine comm ,
  const std::vector<unsigned> & partition ,
  const unsigned iband ,
  const unsigned nband ,
  const unsigned stride ,
  const double evalue ,
  std::vector<unsigned> & prefix ,
  std::vector<unsigned> & coli ,
  std::vector<double>   & coef )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned local_irow = partition[ p_rank ];
  const unsigned local_nrow = partition[ p_rank + 1 ] - local_irow ;

  const unsigned nglobal = partition[ p_size ];
  const unsigned nzrow = 1 + 2 * nband ;
  const unsigned local_alloc = local_nrow * nzrow ;

  prefix.resize( local_nrow + 1 );
  coli.resize( local_alloc );
  coef.resize( local_alloc );

  for ( unsigned i = 0 ; i <= local_nrow ; ++i ) {
    prefix[i] = i * nzrow ;
  }

  std::vector< std::pair<unsigned,double> > row( nzrow );

  for ( unsigned i = 0 ; i < local_nrow ; ++i ) {
    const unsigned irow = local_irow + i ;

    row[0].first  = irow ;
    row[0].second = evalue + ( nzrow - 1 );

    for ( unsigned j = 0 ; j < nband ; ++j ) {
      const unsigned b = iband + j * stride ;
      row[j+1].first  = ( irow + b ) % nglobal ;
      row[j+1].second = -1 ;

      row[j+1+nband].first = ( irow + nglobal - b ) % nglobal ;
      row[j+1+nband].second = -1 ;
    }

    std::sort( row.begin() , row.end() );

    const unsigned k = i * nzrow ;

    for ( unsigned j = 0 ; j < nzrow ; ++j ) {
      coli[ k + j ] = row[j].first ;
      coef[ k + j ] = row[j].second ;
    }
  }
}

struct FillWork {
  double   x_mag ;
  double * x_beg ;
  unsigned x_length ;
};

void task_rand_fill( void * arg , TPI::ThreadPool pool )
{
  const double r_max = RAND_MAX ;
  int p_rank , p_size ; TPI::Pool_rank( pool , p_rank , p_size );

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
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Want to run each test for approximately 1 second.

void timing_txblas1(
  ParallelMachine comm ,
  TPI::ThreadPool pool , 
  const unsigned M ,
  const unsigned Mflop_target )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned m_rem   = M % p_size ;
  const unsigned m_local = M / p_size + ( p_rank < m_rem ? 1 : 0 );

  std::vector<double> x_values( m_local );
  std::vector<double> y_values( m_local );

  if ( x_values.size() != m_local || y_values.size() != m_local ) {
    std::cout << "timing_test_blas1 failed to allocate "
              << m_local << std::endl ;
    return ;
  }

  {
    FillWork data ;
    data.x_mag = 1 ;
    data.x_beg = & x_values[0] ;
    data.x_length = m_local ;
    TPI::Run( pool , & task_rand_fill , & data );
  }

  {
    FillWork data ;
    data.x_mag = 1 ;
    data.x_beg = & y_values[0] ;
    data.x_length = m_local ;
    TPI::Run( pool , & task_rand_fill , & data );
  }

  parallel_machine_barrier( comm );

  //--------------------------------
  // standard dot

  {
    // ddot = ( 1 mult + 1 add ) * M
    const double mflop_cycle = ((double) 2 * M ) / 1.0e6 ;
    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );
    const double mflop_total = mflop_cycle * ncycle ;

    double t = wall_time();

    for ( unsigned i = 0 ; i < ncycle ; ++i ) {
      double s ;
      tddot( pool , & s , m_local , & x_values[0] , & y_values[0] );
      all_reduce( comm , ReduceSum<1>( & s ) );
    }

    double dt = wall_dtime( t );

    all_reduce( comm , ReduceMax<1>( & dt ) );

    if ( p_rank == 0 ) {
      std::cout << "TIMING tddot =   "
                << ( mflop_total / dt ) << " Mflops ( "
                << mflop_cycle << " * " << ncycle << " / " << dt << " )"
                << std::endl ;
      std::cout.flush();
    }
  }

  //--------------------------------
  // double-double accumulation dot

  {
    // xddot = ( 1 mult + 7 add + 2 compare ) * M
    const double mflop_cycle = ((double) 10 * M ) / 1.0e6 ;
    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );
    const double mflop_total = mflop_cycle * ncycle ;

    double t = wall_time();

    for ( unsigned i = 0 ; i < ncycle ; ++i ) {
      Summation z ;
      txddot( pool , z.xdval , m_local , & x_values[0] , & y_values[0] );
      all_reduce( comm , ReduceSum<1>( & z ) );
    }

    double dt = wall_dtime( t );

    all_reduce( comm , ReduceMax<1>( & dt ) );

    if ( p_rank == 0 ) {
      std::cout << "TIMING txddot =  "
                << ( mflop_total / dt ) << " Mflops ( "
                << mflop_cycle << " * " << ncycle << " / " << dt << " )"
                << std::endl ;
      std::cout.flush();
    }
  }

  //--------------------------------
  // standard  y[i] = a * x[i] + b * y[i]

  {
    // daxpby = 2 mult + 1 add
    const double mflop_cycle = ((double) 3 * M ) / 1.0e6 ;
    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );
    const double mflop_total = mflop_cycle * ncycle ;

    double t = wall_time();

    for ( unsigned i = 0 ; i < ncycle ; ++i ) {
      tdaxpby( pool , m_local , 1.0 , & x_values[0] , 1.0 , & y_values[0] );
    }

    double dt = wall_dtime( t );

    all_reduce( comm , ReduceMax<1>( & dt ) );

    if ( p_rank == 0 ) {
      std::cout << "TIMING tdaxpby = "
                << ( mflop_total / dt ) << " Mflops ( "
                << mflop_cycle << " * " << ncycle << " / " << dt << " )"
                << std::endl ;
      std::cout.flush();
    }
  }

  //--------------------------------
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void timing_txblas_cr_mxv( ParallelMachine comm ,
                           TPI::ThreadPool pool ,
                           const unsigned M ,
                           const unsigned nband ,
                           const unsigned stride ,
                           const unsigned Mflop_target ,
                           const double Value )
{
  const unsigned p_rank = parallel_machine_rank( comm );
  const unsigned nz_row = 1 + 2 * nband ;
  const unsigned nz_tot = M * nz_row ;

  std::vector<unsigned> partition ;

  simple_partition( comm , M , partition );

  const unsigned n_local = partition[ p_rank + 1 ] - partition[ p_rank ];
  const unsigned nz_local = n_local * nz_row ;

  const double d_one = 1 ;

  std::vector<double> x_vec( n_local , d_one );
  std::vector<double> y_vec( n_local );

  double * const x = & x_vec[0] ;
  double * const y = & y_vec[0] ;

  std::vector<unsigned> cr_prefix( n_local + 1 );
  std::vector<unsigned> cr_coli( nz_local );
  std::vector<double>   cr_coef( nz_local );

  test_fill_cr_band( comm, partition, 1, nband, stride, Value,
                     cr_prefix , cr_coli , cr_coef );

  CR_Matrix matrix( comm , pool , partition , cr_prefix , cr_coli , cr_coef );

  const double mflop_cycle = ((double) 2 * nz_tot ) / 1.0e6 ;
  const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );
  const double mflop_total = mflop_cycle * ncycle ;

  double t = wall_time();

  for ( unsigned i = 0 ; i < ncycle ; ++i ) {
    matrix.multiply( x , y );
  }

  double dt =  wall_dtime( t );

  all_reduce( comm , ReduceMax<1>( & dt ) );

  if ( p_rank == 0 ) {
    std::cout << "TIMING CR_MXV = " ;
    std::cout << ( mflop_total / dt );
    std::cout << " Mflops ( " ;
    std::cout  << mflop_cycle << " * " << ncycle << " / " << dt << " )";
    std::cout << std::endl ;
  }

  return ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_timing_mxv( phdmesh::ParallelMachine comm ,
                      TPI::ThreadPool pool ,
                      std::istream & is )
{
  unsigned nsize = 1000000 ;
  unsigned nband = 50 ;
  unsigned stride = 1 ;
  unsigned mflop  = 1000 ;
  double   value = 1 ;
  if ( is.good() ) { is >> nsize ; }
  if ( is.good() ) { is >> nband ; }
  if ( is.good() ) { is >> stride ; }
  if ( is.good() ) { is >> mflop ; }

  timing_txblas_cr_mxv(  comm , pool , nsize , nband , stride , mflop , value );
}

void test_timing_blas1( phdmesh::ParallelMachine comm ,
                       TPI::ThreadPool pool ,
                       std::istream & is )
{
  unsigned num = 10000000 ;
  unsigned mflop = 1000 ;
  if ( is.good() ) { is >> num ; }
  if ( is.good() ) { is >> mflop ; }
  timing_txblas1( comm , pool , num , mflop );
}

//----------------------------------------------------------------------


