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

#include <stdio.h>
#include <util/TPI.h>

struct TestTPI {
  int loop_size ;
  int loop_count ;
  int count ;
};

static
void test_tpi_loop( void * arg , TPI_ThreadPool pool )
{
  static const char name[] = "test_tpi_loop" ;

  struct TestTPI * const data = (struct TestTPI *) arg ;

  int nlock  = -1 ;
  int size  = -1 ;
  int rank  = -1 ;
  int result = 0 ;

  if ( ! result && ( result = TPI_Pool_rank( pool , & rank , & size ) ) ) {
    fprintf(stderr,"%s: TPI_Pool_rank = %d\n",name,result);
  }

  if ( ! result && ( result = TPI_Lock_size( pool , & nlock ) ) ) {
    fprintf(stderr,"%s: TPI_Lock_size = %d\n",name,result);
  }

  if ( ! result && nlock < 1 ) {
    fprintf(stderr,"%s: number locks = %d\n",name,nlock);
  }

  if ( ! result && nlock ) {
    const int begin = ( data->loop_size * rank ) / size ;
    const int end   = ( data->loop_size * ( rank + 1 ) ) / size ;

    int count = 0 ;
    for ( int j = 0 ; j < 1001 ; ++j ) {
      if ( j % 2 ) {
        for ( int i = begin ; i < end ; ++i ) { --count ; }
      }
      else {
        for ( int i = begin ; i < end ; ++i ) { ++count ; }
      }
    }

    TPI_Lock( pool , 0 );
    data->count += rank ;
    data->loop_count += count ;
    TPI_Unlock( pool , 0 );
  }

  return ;
}

int test_c_tpi( TPI_ThreadPool pool )
{
  static const char name[] = "test_tpi_loop" ;
  int size = -1 ;
  int result = 0 ;

  if ( ! result && ( result = TPI_Pool_rank( pool , NULL , & size ) ) ) {
    fprintf(stderr,"%s: TPI_Pool_rank = %d\n",name,result);
  }

  if ( ! result ) {
    struct TestTPI data = { 1000000 , 0 , 0 };
    if ( ( result = TPI_Set_lock_size( pool , 1 ) ) ) {
      fprintf(stderr,"%s: TPI_Set_lock_size = %d\n",name,result);
    }
    if ( ( result = TPI_Run( pool , & test_tpi_loop , & data ) ) ) {
      fprintf(stderr,"%s: TPI_Run = %d\n",name,result);
    }
    if ( ! result ) {
      int count = 0 ;
      for ( int i = 0 ; i < size ; ++i ) { count += i ; }

      if ( ( result = data.count != count ) ) {
        fprintf(stderr,"%s: test_tpi_loop : %d != %d\n",name,
                data.count , count );
      }
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

struct TestTPISplit {
  int split_rank ;
  int split_size ;
  int pool_size ;
  int error ;
};


static
void test_tpi_split_run( void * arg , TPI_ThreadPool pool )
{
  struct TestTPISplit * const data = (struct TestTPISplit *) arg ;

  int pool_size ;
  int rank ;
  int result ;

  if ( ( result = TPI_Pool_rank( pool , & rank , & pool_size ) ) ) {
    printf("test_tpi_split failed TPI_Rank = %d\n",result);
  }

  if ( ( result = TPI_Lock( pool , 0 ) ) ) {
    printf("test_tpi_split %d.%d failed TPI_Lock = %d\n",
           data->split_rank , rank , result );
  }

  if ( TPI_Lock( pool , 1 ) ) {
    printf("test_tpi_split %d.%d successfull error check TPI_Lock\n",
           data->split_rank , rank );
  }

  TPI_Unlock( pool , 0 );
}

static
void test_tpi_split( void * arg , TPI_ThreadPool pool )
{
  struct TestTPISplit * const data = (struct TestTPISplit *) arg ;

  int pool_size ;
  int rank ;

  data->error = 0 ;

  if ( TPI_Pool_rank( pool , & rank , & pool_size ) ) {
    data->error = 1 ;
  }
  else if ( pool_size != data->pool_size || 0 != rank ) {
    data->error = 2 ;
  }
  else if ( TPI_Set_lock_size( pool , 1 ) ) {
    data->error = 3 ;
  }
  else if ( TPI_Run( pool , & test_tpi_split_run , arg ) ) {
    data->error = 4 ;
  }
}

int test_c_tpi_split( TPI_ThreadPool pool )
{
  static const char name[] = "test_tpi_split" ;

  int size = -1 ;
  int result = 0 ;

  if ( ! result && ( result = TPI_Pool_rank( pool , NULL , & size ) ) ) {
    fprintf(stderr,"%s: TPI_Pool_rank = %d\n",name,result);
  }

  if ( ! result && 1 < size ) {
    int split_size[ 2 ] = { 2 , size - 2 };
    TPI_parallel_subprogram split_routine[2] =
      { & test_tpi_split , & test_tpi_split };
    struct TestTPISplit data[2] = { { 0 , 2 , 2 , 0 } ,
                                    { 1 , 2 , size - 2 , 0 } };
    void * ptr[2] = { data , data + 1 };

    int nsplit = split_size[1] ? 2 : 1 ;

    result = TPI_Split( pool , nsplit , split_size , split_routine , ptr );

    if      ( data[0].error ) { result = data[0].error ; }
    else if ( data[1].error ) { result = data[1].error ; }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

