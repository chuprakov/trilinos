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

#include <stdexcept>
#include <iostream>
#include <vector>
#include <util/TPI.h>
#include <util/Parallel.hpp>
#include <util/NamedValue.hpp>

extern "C" {

struct TestTPI_Chunk {
  double *  coef ;
  double ** array ;
  unsigned  num_array ;
  unsigned  len_array ;
  unsigned  len_chunk ;
  unsigned  chunk ;
};

static
void test_tpi_chunk_array_locking( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double         * const coef  = data->coef  ;
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
          for ( unsigned j = 0 ; j < num_array ; ++j ) {
            tmp += coef[j] * array[j][i] ;
          }
          array[0][i] = tmp ;
        }
      }
    }
  }

  return ;
}

static
void test_tpi_chunk_array_clean( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double         * const coef  = data->coef  ;
  double * const * const array = data->array  ;
  const unsigned     num_array = data->num_array ;
  const unsigned     len_array = data->len_array ;
  const unsigned     len_chunk = data->len_chunk ;
  const unsigned     num_chunk = len_array / len_chunk +
                               ( len_array % len_chunk ? 1 : 0 );
  int p_rank , p_size ;

  TPI_Pool_rank( pool , & p_rank , & p_size );

  /* array[ num_array ][ len_array ] */

  {
    const unsigned beg_chunk = ( num_chunk * p_rank )         / p_size ;
    const unsigned end_chunk = ( num_chunk * ( p_rank + 1 ) ) / p_size ;

    for ( unsigned chunk = beg_chunk ; chunk < end_chunk ; ++chunk ) {
      const unsigned beg = chunk * len_chunk ;
      const unsigned end = ( beg + len_chunk < len_array )
                           ? beg + len_chunk : len_array ;

      for ( unsigned i = beg ; i < end ; ++i ) {
        double tmp = 0 ;
        for ( unsigned j = 0 ; j < num_array ; ++j ) {
          tmp += coef[j] * array[j][i] ;
        }
        array[0][i] = tmp ;
      }
    }
  }

  return ;
}

static
void test_tpi_chunk_col_locking( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double         * const coef  = data->coef  ;
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
        const unsigned len = end - beg ;
        double * const arr = array[ chunk ];

        for ( unsigned i = 0 ; i < len ; ++i ) {
          double * const a = arr + i ;
          double tmp = 0 ;
          for ( unsigned j = 0 ; j < num_array ; ++j ) {
            tmp += coef[j] * a[j*len_chunk] ;
          }
          a[0] = tmp ;
        }
      }
    }
  }

  return ;
}

static
void test_tpi_chunk_col_clean( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double         * const coef  = data->coef  ;
  double * const * const array = data->array ;
  const unsigned     num_array = data->num_array ;
  const unsigned     len_array = data->len_array ;
  const unsigned     len_chunk = data->len_chunk ;
  const unsigned     num_chunk = len_array / len_chunk +
                               ( len_array % len_chunk ? 1 : 0 );
  int p_rank , p_size ;

  TPI_Pool_rank( pool , & p_rank , & p_size );

  /* array[ num_chunk ][ len_chunk * num_array ] */

  {
    const unsigned beg_chunk = ( num_chunk * p_rank )         / p_size ;
    const unsigned end_chunk = ( num_chunk * ( p_rank + 1 ) ) / p_size ;

    for ( unsigned chunk = beg_chunk ; chunk < end_chunk ; ++chunk ) {
      const unsigned beg = chunk * len_chunk ;
      const unsigned end = ( beg + len_chunk < len_array )
                           ? beg + len_chunk : len_array ;

      const unsigned len = end - beg ;
      double * const arr = array[ chunk ];

      for ( unsigned i = 0 ; i < len ; ++i ) {
        double * const a = arr + i ;
        double tmp = 0 ;
        for ( unsigned j = 0 ; j < num_array ; ++j ) {
          tmp += coef[j] * a[j*len_chunk] ;
        }
        a[0] = tmp ;
      }
    }
  }

  return ;
}

static
void test_tpi_chunk_row_locking( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double         * const coef  = data->coef  ;
  double * const * const array = data->array ;
  const unsigned     num_array = data->num_array ;
  const unsigned     len_array = data->len_array ;
  const unsigned     len_chunk = data->len_chunk ;

  /* array[ num_chunk ][ num_array * len_chunk ] */

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
        const unsigned len = end - beg ;
        const unsigned i_end = len * num_array ;
        double * const arr = array[ chunk ];

        for ( unsigned i = 0 ; i < i_end ; i += num_array ) {
          double * const a = arr + i ;
          double tmp = 0 ;
          for ( unsigned j = 0 ; j < num_array ; ++j ) {
            tmp += coef[j] * a[j] ;
          }
          a[0] = tmp ;
        }
      }
    }
  }

  return ;
}

static
void test_tpi_chunk_row_clean( void * arg , TPI_ThreadPool pool )
{
  struct TestTPI_Chunk * const data = (struct TestTPI_Chunk *) arg ;

  double         * const coef  = data->coef  ;
  double * const * const array = data->array ;
  const unsigned     num_array = data->num_array ;
  const unsigned     len_array = data->len_array ;
  const unsigned     len_chunk = data->len_chunk ;
  const unsigned     num_chunk = len_array / len_chunk +
                               ( len_array % len_chunk ? 1 : 0 );
  int p_rank , p_size ;

  TPI_Pool_rank( pool , & p_rank , & p_size );

  /* array[ num_chunk ][ num_array * len_chunk ] */

  {
    const unsigned beg_chunk = ( num_chunk * p_rank )         / p_size ;
    const unsigned end_chunk = ( num_chunk * ( p_rank + 1 ) ) / p_size ;

    for ( unsigned chunk = beg_chunk ; chunk < end_chunk ; ++chunk ) {
      const unsigned beg = chunk * len_chunk ;
      const unsigned end = ( beg + len_chunk < len_array )
                           ? beg + len_chunk : len_array ;

      const unsigned len = end - beg ;
      double * const arr = array[ chunk ];

      for ( unsigned i = 0 ; i < len ; ++i ) {
        double * const a = arr + i * num_array ;
        double tmp = 0 ;
        for ( unsigned j = 0 ; j < num_array ; ++j ) {
          tmp += coef[j] * a[j] ;
        }
        a[0] = tmp ;
      }
    }
  }

  return ;
}

}

using namespace phdmesh ;

namespace {

void test_tpi_chunk( ParallelMachine , TPI_ThreadPool pool ,
                     const unsigned len_array ,
                     const unsigned len_chunk ,
                     const unsigned locking ,
                     const unsigned Mflop_target )
                   
{
  enum { NUM_TEST = 6 };
  const unsigned test_num_array[ NUM_TEST ] = { 2 , 5 , 10 , 20 , 50 , 100 };

  struct TestTPI_Chunk data = { NULL , NULL , 0 , 0 , 0 , 0 };
  double dt_array[ NUM_TEST ] , mflops_array[ NUM_TEST ];
  double dt_chunk_col[ NUM_TEST ] , mflops_chunk_col[ NUM_TEST ];
  double dt_chunk_row[ NUM_TEST ] , mflops_chunk_row[ NUM_TEST ];
  double t = 0 ;

  data.num_array = 0 ;
  data.len_array = len_array ;
  data.len_chunk = len_chunk ;
  data.coef = (double *) calloc( test_num_array[NUM_TEST-1], sizeof(double) );

  //--------------------------------------------------------------------
  // Array test

  data.array = (double **)
    malloc( sizeof(double*) * test_num_array[NUM_TEST-1] );

  for ( unsigned i = 0 ; i < test_num_array[ NUM_TEST - 1 ] ; ++i ) {
    data.array[i] = (double *) calloc( len_array , sizeof(double) );
    if ( NULL == data.array[i] ) {
      throw std::runtime_error( std::string("Failed to allocate") );
    }
  }

  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {

    const double mflop_cycle =
      ((double)( 2 * test_num_array[ i_test ] * len_array )) / 1.0e6 ;
    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

    data.num_array = test_num_array[ i_test ];

    for ( unsigned repeat = 0 ; repeat < 3 ; ++repeat ) {
      t = wall_time();

      if ( locking ) {
        for ( unsigned i = 0 ; i < ncycle ; ++i ) {
          data.chunk = 0 ;
          TPI_Set_lock_size( pool , 1 );
          TPI_Run( pool , & test_tpi_chunk_array_locking , & data );
        }
      }
      else {
        for ( unsigned i = 0 ; i < ncycle ; ++i ) {
          data.chunk = 0 ;
          TPI_Run( pool , & test_tpi_chunk_array_clean , & data );
        }
      }

      const double dt = wall_dtime( t );

      if ( 0 == repeat || dt < dt_array[ i_test ] ) {
        dt_array[ i_test ] = dt ;
      }
    }
    mflops_array[ i_test ] = mflop_cycle * ncycle / dt_array[ i_test ];
  }

  for ( unsigned i = 0 ; i < test_num_array[ NUM_TEST - 1 ] ; ++i ) {
    free( data.array[i] );
  }
  free( data.array );

  //--------------------------------------------------------------------
  // Chunk tests

  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    const unsigned num_array = test_num_array[ i_test ];
    const unsigned num_chunks = 1 + len_array / len_chunk ;
    const unsigned size_chunk = num_array * len_chunk ;

    const double mflop_cycle =
      ((double)( 2 * num_array * len_array )) / 1.0e6 ;
    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

    data.num_array = num_array ;
    data.array     = (double **) malloc( sizeof(double*) * num_chunks );

    for ( unsigned i = 0 ; i < num_chunks ; ++i ) {
      data.array[i] = (double *) calloc( size_chunk , sizeof(double) );
      if ( NULL == data.array[i] ) {
        throw std::runtime_error( std::string("Failed to allocate") );
      }
    }

    for ( unsigned repeat = 0 ; repeat < 3 ; ++repeat ) {
      t = wall_time();

      if ( locking ) {
        for ( unsigned i = 0 ; i < ncycle ; ++i ) {
          data.chunk = 0 ;
          TPI_Set_lock_size( pool , 1 );
          TPI_Run( pool , & test_tpi_chunk_col_locking , & data );
        }
      }
      else {
        for ( unsigned i = 0 ; i < ncycle ; ++i ) {
          data.chunk = 0 ;
          TPI_Run( pool , & test_tpi_chunk_col_clean , & data );
        }
      }

      const double dt = wall_dtime( t );

      if ( 0 == repeat || dt < dt_chunk_col[ i_test ] ) {
        dt_chunk_col[ i_test ] = dt ;
      }
    }

    mflops_chunk_col[ i_test ] = mflop_cycle * ncycle / dt_chunk_col[ i_test ];

    for ( unsigned repeat = 0 ; repeat < 3 ; ++repeat ) {
      t = wall_time();

      for ( unsigned i = 0 ; i < ncycle ; ++i ) {
        data.chunk = 0 ;
        if ( locking ) {
          TPI_Set_lock_size( pool , 1 );
          TPI_Run( pool , & test_tpi_chunk_row_locking , & data );
        }
        else {
          TPI_Run( pool , & test_tpi_chunk_row_clean , & data );
        }
      }

      const double dt = wall_dtime( t );

      if ( 0 == repeat || dt < dt_chunk_row[ i_test ] ) {
        dt_chunk_row[ i_test ] = dt ;
      }
    }

    mflops_chunk_row[ i_test ] = mflop_cycle * ncycle / dt_chunk_row[ i_test ];

    for ( unsigned i = 0 ; i < num_chunks ; ++i ) {
      free( data.array[i] );
    }
    free( data.array );
  }

  free( data.coef );

  //--------------------------------------------------------------------

  std::cout << std::endl
            << "TIMING MULTIARRAY len_array = " << len_array
            << " , len_chunk = " << len_chunk
            << " , locking = " << locking
            << " , target = " << Mflop_target << " Mflops"
            << std::endl ;
  std::cout << "Arrays =" ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << test_num_array[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;

  std::cout << "FLAT ARRAY test time (sec) = " ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << dt_array[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;
  std::cout << "FLAT ARRAY test Mflops = " ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << mflops_array[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;

  std::cout << "CHUNK COLUMN test time (sec) = " ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << dt_chunk_col[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;
  std::cout << "CHUNK COLUMN test Mflops = " ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << mflops_chunk_col[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;

  std::cout << "CHUNK ROW test time (sec) = " ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << dt_chunk_row[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;
  std::cout << "CHUNK ROW test Mflops = " ;
  for ( unsigned i_test = 0 ; i_test < NUM_TEST ; ++i_test ) {
    std::cout << " " << mflops_chunk_row[ i_test ] ;
  }
  std::cout << std::endl << std::endl ;
  std::cout.flush();

  return ;
}

}

/*--------------------------------------------------------------------*/

void test_tpi_chunk( ParallelMachine machine ,
                     TPI_ThreadPool pool ,
                     std::istream & s )
{
  enum { NUM_CHUNK_SIZE = 7 };
  const unsigned chunk_size[ NUM_CHUNK_SIZE ] =
    { 100 , 200 , 500 , 1000 , 2000 , 5000 , 10000 };

  NamedValue<unsigned> length( "length" );
  NamedValue<unsigned> target( "target" );
  NamedValue<unsigned> locking( "locking" );
  NamedValue< std::vector<unsigned> > chunk( "chunk" );

  length.value = 1e6 ;
  locking.value = 0 ;
  target.value = 200 ;

  {
    NamedValueSet input_values ;

    input_values.insert( target );
    input_values.insert( locking );
    input_values.insert( chunk );

    s >> input_values ;
  }

  if ( chunk.value.empty() ) {
    for ( unsigned i = 0 ; i < NUM_CHUNK_SIZE ; ++i ) {
      test_tpi_chunk( machine , pool ,
                      length.value , chunk_size[i] ,
                      locking.value , target.value );
    }
  }
  else {
    for ( unsigned i = 0 ; i < chunk.value.size() ; ++i ) {
      test_tpi_chunk( machine , pool ,
                      length.value , chunk.value[i] ,
                      locking.value , target.value );
    }
  }

  return ;
}

/*--------------------------------------------------------------------*/

