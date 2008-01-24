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

#include <util/TPI.h>

#include <unistd.h>
#include <errno.h>
#include <pthread.h>

/*--------------------------------------------------------------------*/

struct ThreadPool_Data ;

typedef struct TPI_ThreadPool_Private {
  pthread_mutex_t          m_lock ;
  int                   (* m_check )( TPI_ThreadPool );

  struct ThreadPool_Data * m_pool ;
  void                   * m_buffer ;
} ThreadData ;

typedef struct ThreadPool_Data {
  pthread_mutex_t          m_active ;
  struct ThreadPool_Data * m_parent ;
  ThreadData             * m_thread ;
  int                      m_thread_size ;

  pthread_mutex_t        * m_mutex ;
  int                      m_mutex_size ;
  TPI_parallel_subprogram  m_routine ;
  void                   * m_argument ;
  int                      m_buf_size ;
} ThreadPool ;

/*--------------------------------------------------------------------*/

static int local_thread_pool_check( TPI_ThreadPool local )
{
  return NULL == local ? TPI_ERROR_NULL :
       ( & local_thread_pool_check != local->m_check ? TPI_ERROR_POOL : 0 );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Pool_rank( TPI_ThreadPool local , int * rank , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result ) {
    if ( NULL != rank ) { *rank = local - local->m_pool->m_thread ; }
    if ( NULL != size ) { *size = local->m_pool->m_thread_size ; }
  }
  return result ;
}

int TPI_Buffer( TPI_ThreadPool local , void ** buffer , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result ) {
    if ( NULL != buffer ) { *buffer = local->m_buffer ; }
    if ( NULL != size   ) { *size   = local->m_pool->m_buf_size ; }
  }
  return result ;
}

int TPI_Lock_size( TPI_ThreadPool local , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result && NULL != size ) { *size = local->m_pool->m_mutex_size ; }
  return result ;
}

int TPI_Split_lock_size( TPI_ThreadPool local , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result && NULL != size ) {
    *size = local->m_pool->m_parent
          ? local->m_pool->m_parent->m_mutex_size : 0 ;
  }
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static int local_thread_pool_lock( ThreadPool * pool , int i )
{
  int result = ( NULL != pool && 0 <= i && i < pool->m_mutex_size )
               ? 0 : TPI_ERROR_SIZE ;
  if ( ! result ) {
    pthread_mutex_t * m = pool->m_mutex + i ;
#if 0
    while ( EBUSY == ( result = pthread_mutex_trylock(m) ) ); 
    if ( result ) { result = TPI_ERROR_LOCK ; }
#else
    result = pthread_mutex_lock(m) ? TPI_ERROR_LOCK : 0 ;
#endif
  }
  return result ;
}

static int local_thread_pool_trylock( ThreadPool * pool , int i )
{
  int result = ( NULL != pool && 0 <= i && i < pool->m_mutex_size )
               ? 0 : TPI_ERROR_SIZE ;
  if ( ! result ) {
    result = pthread_mutex_trylock( pool->m_mutex + i );
    if ( EBUSY == result ) { result = TPI_LOCK_BUSY ; }
    else if ( result )     { result = TPI_ERROR_LOCK ; }
  }
  return result ;
}

static int local_thread_pool_unlock( ThreadPool * pool , int i )
{
  int result = ( NULL != pool && 0 <= i && i < pool->m_mutex_size )
               ? 0 : TPI_ERROR_SIZE ;
  if ( ! result ) {
    result = pthread_mutex_unlock( pool->m_mutex + i ) ? TPI_ERROR_LOCK : 0 ;
  }
  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Lock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) { result = local_thread_pool_lock( local->m_pool , i ); }
  return result ;
}

int TPI_Split_lock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_lock( local->m_pool->m_parent , i );
  }
  return result ;
}

int TPI_Trylock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_trylock( local->m_pool , i );
  }
  return result ;
}

int TPI_Split_trylock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_trylock( local->m_pool->m_parent , i );
  }
  return result ;
}

int TPI_Unlock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) { result = local_thread_pool_unlock( local->m_pool , i ); }
  return result ;
}

int TPI_Split_unlock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_unlock( local->m_pool->m_parent , i );
  }
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static int local_thread_pool_run_routine( ThreadData * const td )
{
  int working = 1 ;

  const ThreadPool * const p = td->m_pool ;

  if ( NULL != p ) {
    const unsigned size = p->m_buf_size ;

    if ( ( working = ( NULL != p->m_routine ) ) ) {
      char buffer[ size ];
      void * const old = td->m_buffer ;
      td->m_buffer = size ? buffer : NULL ;
      (* p->m_routine)( p->m_argument, td );
      td->m_buffer = old ;
    }
    td->m_pool = NULL ; /* Have completed */
  }

  return working ;
}

/*--------------------------------------------------------------------*/
/* Driver: call the routine, check for thread termination */

static void * local_thread_pool_driver( void * arg )
{
  ThreadData * const td = (ThreadData*)( arg );

  for ( int working = 1 ; working ; ) {
#if 0
    /* spin lock */
    if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
      working = local_thread_pool_run_routine( td );
      pthread_mutex_unlock( & td->m_lock );
    }
#else
    /* hard lock */
    pthread_mutex_lock( & td->m_lock );
    working = local_thread_pool_run_routine( td );
    pthread_mutex_unlock( & td->m_lock );
#endif
  }

  return NULL ;
}

/*--------------------------------------------------------------------*/

static int local_thread_pool_run( ThreadPool * const pool )
{
  const unsigned number_thread = pool->m_thread_size ;
  const unsigned number_locks  = pool->m_mutex_size ;

  int result = 0 ;

  {
    pthread_mutex_t locks[ number_locks ];

    unsigned nlocks = 0 ;

    while ( nlocks < number_locks && ! result ) {
      if ( pthread_mutex_init( locks + nlocks , NULL ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        ++nlocks ;
      }
    }

    if ( ! result ) {
      int running[ number_thread ];

      pool->m_mutex = locks ;

      for ( unsigned i = 1 ; i < number_thread ; ++i ) {
        ThreadData * const td = pool->m_thread + i ;
        if ( ( running[i] = NULL != td->m_pool ) ) {
          pthread_mutex_unlock( & td->m_lock );
        }
      }

      local_thread_pool_run_routine( pool->m_thread );

      for ( int waiting = 1 ; waiting ; ) {
        waiting = 0 ;
        for ( unsigned i = 1 ; i < number_thread ; ++i ) {
          if ( running[i] ) {
            ThreadData * const td = pool->m_thread + i ;
            if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
              local_thread_pool_run_routine( td );
              running[i] = 0 ;
            }
            else {
              waiting = 1 ;
            }
          }
        }
      }
    }

    while ( nlocks-- ) { pthread_mutex_destroy( locks + nlocks ); }
  }

  for ( unsigned i = 0 ; i < number_thread ; ++i ) {
    pool->m_thread[i].m_pool = NULL ;
  }

  pool->m_mutex      = NULL ;
  pool->m_mutex_size = 0 ;
  pool->m_routine    = NULL ;
  pool->m_argument   = NULL ;
  pool->m_buf_size   = 0 ;

  return result ;
}

/*--------------------------------------------------------------------*/
/*  Acquire the active lock for the pool.
 *  This thread must be the first thread in the pool.
 */

static int local_thread_pool_acquire_active_lock( TPI_ThreadPool local )
{
  int result = local_thread_pool_check( local );

  if ( ! result ) {
    ThreadPool * const pool = local->m_pool ;
    if ( local != pool->m_thread || pthread_mutex_trylock(& pool->m_active) ) {
      result = TPI_ERROR_ACTIVE ;
    }
  }
  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Set_lock_size( TPI_ThreadPool local , int size )
{
  int result = local_thread_pool_acquire_active_lock( local );
  if ( ! result ) {
    local->m_pool->m_mutex_size = size ;
    pthread_mutex_unlock( & local->m_pool->m_active );
  }
  return result ;
}

int TPI_Set_buffer_size( TPI_ThreadPool local , int size )
{
  int result = local_thread_pool_acquire_active_lock( local );
  if ( ! result ) {
    local->m_pool->m_buf_size = size ;
    pthread_mutex_unlock( & local->m_pool->m_active );
  }
  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Run( TPI_ThreadPool local ,
             TPI_parallel_subprogram routine ,
             void * routine_data )
{
  int result = local_thread_pool_acquire_active_lock( local );

  if ( ! result ) { /* Now have the active lock for the whole pool */

    ThreadPool * const pool = local->m_pool ;

    pool->m_routine  = routine ;
    pool->m_argument = routine_data ;

    const unsigned number_thread = pool->m_thread_size ;

    for ( unsigned i = 1 ; i < number_thread ; ++i ) {
      pool->m_thread[i].m_pool = pool ;
    }

    result = local_thread_pool_run( pool );

    local->m_pool = pool ;

    pthread_mutex_unlock( & pool->m_active );
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Split( TPI_ThreadPool local ,
               const int number ,
               const int sizes[] ,
               TPI_parallel_subprogram routine[] ,
               void * routine_data[] )
{
  int result = local_thread_pool_acquire_active_lock( local );

  if ( ! result ) {

    ThreadPool * const parent = local->m_pool ;

    {
      int child_size = 0 ;

      for ( int i = 0 ; i < number && ! result ; ++i ) {
        if      ( NULL == routine[i] ) { result = TPI_ERROR_NULL ; }
        else if ( 0    >= sizes[i] )   { result = TPI_ERROR_SIZE ; }
        child_size += sizes[i] ;
      }

      if ( ! result && child_size != parent->m_thread_size ) {
        result = TPI_ERROR_SIZE ;
      }
    }

    if ( ! result ) {

      ThreadPool children[ number ];
      int child = 0 ;

      for ( int offset = 0 ; child < number && ! result ; ) {

        ThreadData * const td = parent->m_thread + offset ;

        offset += sizes[child] ;

        ThreadPool * const pool = children + child ;

        td->m_pool = pool ;

        pool->m_parent      = parent ;
        pool->m_thread      = td ;
        pool->m_thread_size = sizes[child] ;
        pool->m_mutex       = NULL ;
        pool->m_mutex_size  = 0 ;
        pool->m_routine     = routine[child] ;
        pool->m_argument    = routine_data[child] ;
        pool->m_buf_size    = 0 ;

        if ( pthread_mutex_init( & pool->m_active , NULL ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
        else {
          ++child ;
        }
      }

      if ( ! result ) { result = local_thread_pool_run( parent ); }

      while ( child-- ) {
        ThreadPool * const pool = children + child ;
        pthread_mutex_destroy( & pool->m_active );
        pool->m_thread->m_pool = NULL ;
      }

      local->m_pool = parent ;
    }

    pthread_mutex_unlock( & parent->m_active );
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void local_thread_pool_start_task( void * arg , TPI_ThreadPool local )
{
  if ( NULL  != arg &&
       NULL  != local &&
       NULL  != local->m_pool &&
       local == local->m_pool->m_thread + *( (int*) arg ) ) {
    local->m_check = & local_thread_pool_check ;
  }
}

static void local_thread_pool_terminate_threads( ThreadPool * pool )
{
  const int n = pool->m_thread_size ;

  pool->m_thread_size = 0 ;
  pool->m_mutex_size = 0 ;
  pool->m_routine  = NULL ;
  pool->m_argument = NULL ;
  pool->m_buf_size = 0 ;

  for ( int i = 1 ; i < n ; ++i ) {
    ThreadData * const td = pool->m_thread + i ;
    td->m_pool = pool ;
    pthread_mutex_unlock( & td->m_lock );

    for ( int waiting = 1 ; waiting ; ) {
      if ( ! pthread_mutex_lock( & td->m_lock ) ) {
        if ( NULL == td->m_pool ) { waiting = 0 ; }
        pthread_mutex_unlock( & td->m_lock );
      }
    }
    pthread_mutex_destroy( & td->m_lock );
  }
}

static int local_thread_pool_create_threads( int n , ThreadPool * pool )
{
  int result = pool->m_thread_size ? TPI_ERROR_ACTIVE : 0 ;
  int n_thread = 0 ;
  pthread_attr_t thread_attr ;

  if ( ! result && pthread_attr_init( & thread_attr ) ) {
    result = TPI_ERROR_INTERNAL ;
  }

  if ( ! result ) {
    /* pthread_setconcurrency( n ); */

    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

/* GNU runtime rejects this attribute:
    pthread_attr_setschedpolicy( & thread_attr, SCHED_FIFO );
*/

    pool->m_routine  = & local_thread_pool_start_task ;
    pool->m_argument = & n_thread ;

    pool->m_thread->m_check  = & local_thread_pool_check ;
    pool->m_thread->m_pool   = pool ;
    pool->m_thread->m_buffer = NULL ;

    for ( n_thread = 1 ; n_thread < n && ! result ; ) {
      pthread_t pt ;

      ThreadData * const td = pool->m_thread + n_thread ;

      td->m_check  = NULL ;
      td->m_pool   = pool ;
      td->m_buffer = NULL ;

      if ( pthread_mutex_init( & td->m_lock , NULL ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else if ( pthread_create( & pt ,
                                & thread_attr ,
                                & local_thread_pool_driver ,
                                (void *) td ) ) {
        pthread_mutex_destroy( & td->m_lock );
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        for ( int waiting = 1 ; waiting && ! result ; ) {
          if ( ! pthread_mutex_lock( & td->m_lock ) ) {
            if ( NULL != td->m_pool ) {
              pthread_mutex_unlock( & td->m_lock ); /* has not run */
            }
            else if ( & local_thread_pool_check != td->m_check ) {
              result = TPI_ERROR_INTERNAL ;
            }
            else {
              waiting = 0 ;
              ++n_thread ;
            }
          }
        }
      }
    }
    pthread_attr_destroy( & thread_attr );

    pool->m_thread_size = n_thread ;
    pool->m_routine  = NULL ;
    pool->m_argument = NULL ;

    if ( result ) { local_thread_pool_terminate_threads( pool ); }
  }

  return result ;
}

static int local_thread_pool_root( int n , TPI_ThreadPool * local )
{
  enum { MAX_THREADS = 128 };

  static ThreadData threads[ MAX_THREADS ];

  static ThreadPool root = {
      /* m_active      */  PTHREAD_MUTEX_INITIALIZER ,
      /* m_parent      */  NULL ,
      /* m_thread      */  threads ,
      /* m_thread_size */  0 ,
      /* m_mutex       */  NULL ,
      /* m_mutex_size  */  0 ,
      /* m_routine     */  NULL ,
      /* m_argument    */  NULL ,
      /* m_buf_size    */  0 };

  int result = 0 < n && n < MAX_THREADS ? 0 : TPI_ERROR_SIZE ;

  if ( ! result && pthread_mutex_trylock( & root.m_active ) ) {
    result = TPI_ERROR_ACTIVE ;
  }

  if ( ! result ) {
    if ( local ) {
      result = local_thread_pool_create_threads( n , & root );
      *local = result ? NULL : threads ;
    }
    else {
      local_thread_pool_terminate_threads( & root );
    }
    pthread_mutex_unlock( & root.m_active );
  }

  return result ;
}

int TPI_Init( int n , TPI_ThreadPool * local )
{
  int result = NULL != local ? 0 : TPI_ERROR_NULL ;
  if ( ! result ) { result = local_thread_pool_root( n , local ); }
  return result ;
}

int TPI_Finalize()
{
  return local_thread_pool_root( 1 , NULL );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

