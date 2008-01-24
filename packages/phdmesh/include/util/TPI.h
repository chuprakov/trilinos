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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 */

#ifndef ThreadPoolInterface_h
#define ThreadPoolInterface_h

#if defined( __cplusplus )
extern "C" {
#endif

/*--------------------------------------------------------------------*/

struct TPI_ThreadPool_Private ; 

typedef struct TPI_ThreadPool_Private * TPI_ThreadPool ;

/*--------------------------------------------------------------------*/
/* All functions return zero for success. */

#define TPI_LOCK_BUSY      ((int)  1  /* trylock or unlock failed */ )
#define TPI_ERROR_NULL     ((int) -1  /* NULL input */ )
#define TPI_ERROR_SIZE     ((int) -2  /* BAD input: size or index */ )
#define TPI_ERROR_POOL     ((int) -3  /* BAD input: the pool is invalid */ )
#define TPI_ERROR_LOCK     ((int) -4  /* BAD lock or unlock */ )
#define TPI_ERROR_ACTIVE   ((int) -5  /* BAD input: the pool is active  */ )
#define TPI_ERROR_INTERNAL ((int) -6  /* internal resource error */ )

/*--------------------------------------------------------------------*/
/** Initialize the root thread pool to the specified size. */
int TPI_Init( int , TPI_ThreadPool * );

/** Finalize (shut down) all threads and thread pools. */
int TPI_Finalize();

int TPI_Pool_rank( TPI_ThreadPool , int * /* rank */ , int * /* size */ );

int TPI_Buffer( TPI_ThreadPool , void ** /* buffer */ , int * /* byte_size */ );

int TPI_Lock_size( TPI_ThreadPool , int * );

/* Set the lock size for the next call to TPI_Run. */
int TPI_Set_lock_size( TPI_ThreadPool , int );

/* Set the stack buffer size for the next call to TPI_Run. */
int TPI_Set_buffer_size( TPI_ThreadPool , int );

/*--------------------------------------------------------------------*/
/**  A thread-pool parallel subprogram and its shared data
 *   running within a thread-pool of 'size' threads where
 *   this is the 'rank' thread.
 */
typedef void (*TPI_parallel_subprogram)( void * shared_data ,
                                         TPI_ThreadPool pool );

/** Run a thread-pool parallel subprogram.
 *  Each thread in the pool will call the subprogram as:
 *
 *    (*subprogram)( shared_data , pool )
 *
 *  Nested calls to this routine are illegal.
 */
int TPI_Run( TPI_ThreadPool ,
             TPI_parallel_subprogram ,
             void * /* shared data */ );

/** Blocks until lock # is obtained */
int TPI_Lock( TPI_ThreadPool , int );

/** Tries to lock #, returns 0 if successful */
int TPI_Trylock( TPI_ThreadPool , int );

/** Unlocks lock # */
int TPI_Unlock( TPI_ThreadPool , int );

/*--------------------------------------------------------------------*/
/** Split a thread pool and call the given thread-pool parallel
 *  subprogram on the root thread of each thread-pool with the subset
 *  thread-pool object.
 *
 *    subprogram[i]( data[i] , pool , pool_rank );
 *
 *  These new 'main' subprograms may now call TPI_Run within their child pool.
 */
int TPI_Split( TPI_ThreadPool ,
               const int                  /* Number of child pools     */ ,
               const int []               /* array of child pool sizes */ ,
               TPI_parallel_subprogram [] /* array of main functions  */ ,
               void * []                  /* array of main functions' data */);

int TPI_Split_lock_size( TPI_ThreadPool , int * );
int TPI_Split_lock(      TPI_ThreadPool , int );
int TPI_Split_trylock(   TPI_ThreadPool , int );
int TPI_Split_unlock(    TPI_ThreadPool , int );

/*--------------------------------------------------------------------*/

#if defined( __cplusplus )
}
#endif

#endif

