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

#include <iostream>
#include <vector>
#include <util/TPI.h>
#include <util/Parallel.hpp>
#include <util/NamedValue.hpp>

extern "C" {

struct TestTPI_Chunk {
  double ** array ;
  unsigned  num_array ;
  unsigned  len_array ;
  unsigned  len_chunk ;
  unsigned  chunk ;
};

static
void test_tpi_chunk_array( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double * const * const array = data->array  ;
  const unsigned     num_array = data->num_array ;
  const unsigned     len_array = data->len_array ;
  const unsigned     len_chunk = data->len_chunk ;

  /* array[ num_array ][ len_array ] */

  for ( int work = 1 ; work ; ) {
    unsigned chunk ;

    TPI_Lock( pool , 0 );
    chunk = data->chunk ; ++(data->chunk);
    TPI_Unlock( pool , 0 );

    {
      const unsigned beg = chunk * len_chunk ;
      const unsigned end = ( beg + len_chunk < len_array )
                           ? beg + len_chunk : len_array ;

      if ( ( work = beg < end ) ) {
        for ( unsigned i = beg ; i < end ; ++i ) {
          double tmp = 0 ;
          for ( unsigned j = 1 ; j < num_array ; ++j ) {
            tmp += array[j][i] ;
          }
          array[0][i] += tmp ;
        }
      }
    }
  }

  return ;
}

static
void test_tpi_chunk_chunk( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double * const * const array = data->array ;
  const unsigned     num_array = data->num_array ;
  const unsigned     len_array = data->len_array ;
  const unsigned     len_chunk = data->len_chunk ;

  /* array[ num_chunk ][ len_chunk * num_array ] */

  for ( int work = 1 ; work ; ) {
    unsigned chunk ;

    TPI_Lock( pool , 0 );
    chunk = data->chunk ; ++(data->chunk);
    TPI_Unlock( pool , 0 );

    {
      const unsigned beg = chunk * len_chunk ;
      const unsigned end = ( beg + len_chunk < len_array )
                           ? beg + len_chunk : len_array ;

      if ( ( work = beg < end ) ) {
        const unsigned num   = end - beg ;
        const unsigned j_end = num * num_array ;
        double * const arr = array[ chunk ];

        for ( unsigned i = 0 ; i < num ; ++i ) {
          double * const a = arr + i ;
          double tmp = 0 ;
          for ( unsigned j = num ; j < j_end ; j += num ) {
            tmp += a[j] ;
          }
          a[0] += tmp ;
        }
      }
    }
  }

  return ;
}

}

using namespace phdmesh ;

void test_tpi_chunk( ParallelMachine , TPI_ThreadPool pool , std::istream & s )
{
  struct TestTPI_Chunk data = { NULL , 0 , 0 , 0 , 0 };
  unsigned nloop = 1 ;
  double t = 0 ;
  double dt_array = 0 ;
  double dt_chunk = 0 ;

  {
    NamedValue<unsigned &> num_array( "num_array", data.num_array );
    NamedValue<unsigned &> len_array( "len_array", data.len_array );
    NamedValue<unsigned &> len_chunk( "len_chunk", data.len_chunk );
    NamedValue<unsigned &> num_loop(  "loop" , nloop );

    NamedValueSet input_values ;

    input_values.insert( num_array );
    input_values.insert( len_array );
    input_values.insert( len_chunk );
    input_values.insert( num_loop );

    s >> input_values ;
  }

  if ( s.good() ) {
    { // Chunks
      const unsigned num_chunks = 1 + data.len_array / data.len_chunk ;

      data.array = (double **) malloc( sizeof(double*) * num_chunks );

      for ( unsigned i = 0 ; i < num_chunks ; ++i ) {
        const unsigned beg = i * data.len_chunk ;
        const unsigned end = ( beg + data.len_chunk < data.len_array )
                             ? beg + data.len_chunk : data.len_array ;
        const unsigned num = end - beg ;

        data.array[i] = (double *)
          calloc( num * data.num_array , sizeof(double) );
      }

      data.chunk = 0 ;
      TPI_Set_lock_size( pool , 1 );
      TPI_Run( pool , & test_tpi_chunk_chunk , & data );

      t = wall_time();

      for ( unsigned i = 0 ; i < nloop ; ++i ) {
        data.chunk = 0 ;
        TPI_Set_lock_size( pool , 1 );
        TPI_Run( pool , & test_tpi_chunk_chunk , & data );
      }

      dt_chunk = wall_dtime( t ) / nloop ;

      for ( unsigned i = 0 ; i < num_chunks ; ++i ) {
        free( data.array[i] );
      }
      free( data.array );
    }

    { // Array
      data.array = (double **) malloc( sizeof(double*) * data.num_array );

      for ( unsigned i = 0 ; i < data.num_array ; ++i ) {
        data.array[i] = (double *) calloc( data.len_array , sizeof(double) );
      }

      data.chunk = 0 ;
      TPI_Set_lock_size( pool , 1 );
      TPI_Run( pool , & test_tpi_chunk_array , & data );

      t = wall_time();

      for ( unsigned i = 0 ; i < nloop ; ++i ) {
        data.chunk = 0 ;
        TPI_Set_lock_size( pool , 1 );
        TPI_Run( pool , & test_tpi_chunk_array , & data );
      }

      dt_array = wall_dtime( t ) / nloop ;

      for ( unsigned i = 0 ; i < data.num_array ; ++i ) {
        free( data.array[i] );
      }
      free( data.array );
    }
  }

  std::cout << "TEST_TPI_CHUNK, relative speed = "
            << ( dt_chunk / dt_array )
            << std::endl
            << "    ARRAY DT = " << dt_array
            << " , DT/N = " << ( dt_array / data.num_array )
            << std::endl
            << "    CHUNK DT = " << dt_chunk
            << " , DT/N = " << ( dt_chunk / data.num_array )
            << std::endl ;

  return ;
}

/*--------------------------------------------------------------------*/

