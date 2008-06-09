/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
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
#include <TPI.h>


/*------------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void test_tpi_noop( void * , TPI_ThreadPool );

int test_c_tpi_noop( int num_test , int * num_thread )
{
  const unsigned n = 1e5 ;
  const unsigned n_trial = 7 ;
  int itest ;

  fprintf(stdout,"\n\"TPI_Run(noop) test\"\n");
  fprintf(stdout,"\"NUMBER OF TRIALS\" , %u\n", n_trial );
  fprintf(stdout,
          "\"# Thread\" , \"Min microsec\" , \"Mean microsec\" , \"Max microsec\"\n");

  for ( itest = 0 ; itest < num_test ; ++itest ) {
    const unsigned num = num_thread[ itest ];
    const unsigned n_loop = n / num ;

    unsigned i , j ;
    double dt_min = 0 , dt_max = 0 , dt_mean = 0 ;
    double dt ;

    /* Time many tries and trials, get the min and max time */

    for ( j = 0 ; j < n_trial ; ++j ) {
  
      TPI_Init( num );

      dt = TPI_Walltime();
      for ( i = 0 ; i < n_loop ; ++i ) { TPI_Run( & test_tpi_noop, NULL, 0 ); }
      dt = TPI_Walltime() - dt ;

      dt_mean += dt ;

      if ( ! j ) {
        dt_min = dt_max = dt ;
      }

      if ( dt < dt_min ) { dt_min = dt ; }
      if ( dt > dt_max ) { dt_max = dt ; }

      TPI_Finalize();
    }

    dt_min  *= 1.0e6 / n_loop ;
    dt_mean *= 1.0e6 / ( n_loop * n_trial );
    dt_max  *= 1.0e6 / n_loop ;

    fprintf(stdout, "%u , %g , %g , %g\n", num, dt_min, dt_mean, dt_max );

    fflush(stdout);
  }

  return 0 ;
}

static void test_tpi_noop( void * arg , TPI_ThreadPool pool )
{ while ( arg && pool ) { arg = *((void**)arg) ; } return ; }

/*------------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

struct TestTPI {
  int total ;
  int count ;
};

static void test_tpi_loop( void * arg , TPI_ThreadPool pool )
{
  static const char name[] = "test_tpi_loop" ;

  struct TestTPI * const data = (struct TestTPI *) arg ;

  int size  = -1 ;
  int rank  = -1 ;
  int result = 0 ;
  int begin  = 0 ;
  int number = 0 ;

  if ( ( result = TPI_Rank( pool , & rank , & size ) ) ) {
    fprintf(stderr,"\n%s: TPI_Pool_rank = %d\n",name,result);
  }

  if ( ( result = TPI_Partition(rank,size,data->total, & begin, & number) ) ) {
    fprintf(stderr,"\n%s: TPI_Partition = %d\n",name,result);
  }
  else {
    int count = 0 ;
    int i , j ;
    for ( j = 0 ; j < 101 ; ++j ) {
      if ( j % 2 ) {
        for ( i = 0 ; i < number ; ++i ) { --count ; }
      }
      else {
        for ( i = 0 ; i < number ; ++i ) { ++count ; }
      }
    }

    TPI_Lock( pool , 0 );
    data->count += count ;
    TPI_Unlock( pool , 0 );
  }

  return ;
}

int test_c_tpi_single( int size )
{
  static const char name[] = "test_c_tpi_single" ;
  int result = 0 ;

  struct TestTPI data = { 10000 /* 1000000 */ , 0 };

  TPI_Init( size );

  fprintf(stdout,"\"%s[%d] starting...",name,size);
  fflush(stdout);

  /*--------------------------------*/

  {
    int n ;
    for ( n = 1 ; n < 64 ; ++n ) {
      data.count = 0 ;
      if ( ( result = TPI_Set_lock_size( 1 ) ) ) {
        fprintf(stderr,"\n%s: TPI_Set_lock_size = %d\n",name,result);
      }
      if ( ( result = TPI_Run( & test_tpi_loop , & data , 0 ) ) ) {
        fprintf(stderr,"\n%s: TPI_Run = %d\n",name,result);
      }
      else {
        if ( ( result = data.count != data.total ) ) {
          fprintf(stderr,"\n%s: test_tpi_loop : %d != %d\n",name,
                  data.count , data.total );
        }
      }
    }
  }

  /*--------------------------------*/

  fprintf(stdout,"completed successfully\"\n");
  fflush(stdout);

  TPI_Finalize();

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

