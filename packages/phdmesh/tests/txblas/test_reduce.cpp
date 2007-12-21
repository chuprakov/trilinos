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
#include <util/TPI.h>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <txblas/Reduction.hpp>

using namespace phdmesh ;

void test_reduce( ParallelMachine comm , TPI_ThreadPool pool ,
                  std::istream & is )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  unsigned M = 1000000 ;
  double mag = 1e10 ;

  if ( is.good() ) { is >> M ; }
  if ( is.good() ) { is >> mag ; }

  //--------------------------------------------------------------------
  {
    srand( p_rank ); // each processor do a different series

    const double r_max = RAND_MAX ;

    std::vector<double> values( M );

    for ( std::vector<double>::iterator
          i = values.begin() ; i != values.end() ; ++i ) {
      const double r = rand() ;
      const double x = 2.0 * ( ( r / r_max ) - 0.5 );
      *i = x * mag ;
    }

    std::vector<double> values_neg( values );

    for ( std::vector<double>::iterator
          i = values_neg.begin() ; i != values_neg.end() ; ++i ) {
      *i = - *i ;
    }

    const double a0 = 1 ;
    double a = a0 ;
    double z[4] = { 0 , 0 , 0 , 0 };
    xdsum_add_value( z , a0 );

    for ( std::vector<double>::iterator
          i = values.begin() ; i != values.end() ; ++i ) {
      a += *i ;
    }

    txdsum_add_array( pool , z , M , & values[0] );

    for ( std::vector<double>::iterator
          i = values_neg.begin() ; i != values_neg.end() ; ++i ) {
      a += *i ;
    }

    txdsum_add_array( pool , z , M , & values_neg[0] );

    xdsum_get_value( z , z );

    const int local_flag = z[0] < a0 || a0 < z[0] ;
    int flag = local_flag ;

    all_reduce( comm , ReduceMax<1>( & flag ) );

    if ( flag ) {
      if ( local_flag ) {
        std::cout << "P" << p_rank << " SUMMATION ERROR = "
                  << ( z[0] - a0 ) << std::endl ;
      }
      throw std::runtime_error( std::string( "FAIL SUMMATION TEST 1" ) );
    }

    if ( p_rank == 0 ) {
      std::cout << "PASS SUMMATION TEST: P0 has "
                << z[0] << " vs. " << a << std::endl ;
    }
  }
  //--------------------------------------------------------------------

  double dval[5] = { 0 , 1 , 2 , 3 , 0 };
  long   ival[3] = { 0 , 0 , 0 };
  Summation aval ;
  unsigned flag  = 0 ;

  ival[1] = p_rank + 1 ;
  ival[2] = - ((long)( p_rank + 1 ));
  dval[4] = p_rank + 1 ;
  aval += dval[1] ;

  all_reduce( comm , ReduceSum<5>( dval ) &
                     ReduceMax<3>( ival ) &
                     ReduceSum<1>( & aval ) &
                     ReduceBitOr<1>( & flag ) );

  const double sum_ranks = ( p_size * ( p_size + 1 ) ) / 2 ;

  if ( (unsigned) dval[0] != 0 ||
       (unsigned) dval[1] != p_size ||
       (unsigned) dval[2] != p_size * 2 ||
       (unsigned) dval[3] != p_size * 3 ||
       (unsigned) dval[4] != (unsigned) sum_ranks ||
       ival[0] != 0 ||
       ival[1] != ( (long) p_size ) ||
       ival[2] != -1 ||
       (unsigned) ( aval.value() ) != p_size ||
       flag != 0 ) {
    std::cout << "P" << p_rank << " "
              << dval[0] << " " << dval[1] << " "
              << dval[2] << " " << dval[3] << " "
              << dval[4] << " , "
              << ival[0] << " " << ival[1] << " "
              << ival[2] << " , "
              << aval.value() << " , "
              << flag << std::endl ;
    throw std::runtime_error( std::string( "FAIL REDUCE TEST 1" ) );
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS REDUCE TEST 1 for NP = " << p_size << " : "
              << dval[0] << " " << dval[1] << " "
              << dval[2] << " " << dval[3] << " "
              << dval[4] << " , "
              << ival[0] << " " << ival[1] << " "
              << ival[2] << " , " << flag << std::endl ;
    std::cout.flush();
  }
}

