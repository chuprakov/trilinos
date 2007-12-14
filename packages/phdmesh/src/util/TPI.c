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

#include <phdmesh_config.h>
#include <util/TPI.h>

enum { MAX_THREADS = 1024 };
enum { MAX_LOCK = 15 };

#if defined( PHDMESH_HAS_PTHREADS )

#include <pthread.h>
#include <unistd.h>
#include <errno.h>

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef struct ThreadDataStruct {
  pthread_t               m_thread ;
  pthread_mutex_t         m_lock ;
  TPI_ThreadPool          m_pool ;
  int                     m_size ;
  int                     m_rank ;
  TPI_parallel_subprogram m_routine ;
  void                  * m_argument ;
} ThreadData ;

typedef struct TPI_ThreadPool_Private {
  ThreadData    * m_begin ;
  ThreadData    * m_end ;
  int             m_mutex_size ;
  pthread_mutex_t m_mutex[ MAX_LOCK + 1 ];
} ThreadPool ;

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static TPI_ThreadPool local_thread_pool_root()
{
  static ThreadPool root_pool ;
  return & root_pool ;
}

int TPI_Pool_size( TPI_ThreadPool pool , int * size )
{
  const TPI_ThreadPool root = local_thread_pool_root();

  int result = pool != NULL ? 0 : TPI_ERROR_NULL ;

  if ( ! result ) {

    if ( pool->m_begin < root->m_begin || root->m_end < pool->m_begin ||
         pool->m_end   < root->m_begin || root->m_end < pool->m_end ) {
      result = TPI_ERROR_POOL ;
    }
    else if ( size ) {
      *size = pool->m_end - pool->m_begin ;
    }
  }

  return result ;
}

static int local_thread_pool_lock_create( TPI_ThreadPool pool )
{
  return pthread_mutex_init( pool->m_mutex + MAX_LOCK , NULL )
         ? TPI_ERROR_INTERNAL : 0 ;
}

static void local_thread_pool_lock_destroy( TPI_ThreadPool pool )
{
  pthread_mutex_destroy( pool->m_mutex + MAX_LOCK );
}

static int local_thread_pool_lock( TPI_ThreadPool pool )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result && pthread_mutex_trylock( pool->m_mutex + MAX_LOCK ) ) {
    result = TPI_ERROR_ACTIVE ;
  }

  return result ;
}

static void local_thread_pool_unlock( TPI_ThreadPool pool )
{
  pthread_mutex_unlock( pool->m_mutex + MAX_LOCK );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Lock_allocation( TPI_ThreadPool pool , int number )
{
  int result = local_thread_pool_lock( pool );

  if ( ! result ) {

    local_thread_pool_unlock( pool );

    if ( MAX_LOCK < number ) {
      result = TPI_ERROR_SIZE ;
    }
    else if ( number < 0 ) {
      for ( int i = 0 ; i < pool->m_mutex_size ; ++i ) {
        pthread_mutex_destroy( pool->m_mutex + i );
      }
      pool->m_mutex_size = 0 ;
    }
    else if ( pool->m_mutex_size < number ) {
      for ( int i = pool->m_mutex_size ; i < number ; ++i ) {
        if ( pthread_mutex_init( pool->m_mutex + i , NULL ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
      }
      pool->m_mutex_size = number ;
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Lock_size( TPI_ThreadPool pool , int * size )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result && size ) { *size = pool->m_mutex_size ; }

  return result ;
}

int TPI_Lock( TPI_ThreadPool pool , int i )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result ) {
    result = 0 <= i && i < pool->m_mutex_size ? 0 : TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    pthread_mutex_t * const m = pool->m_mutex + i ;
#if 1
    for ( result = EBUSY ; EBUSY == ( result = pthread_mutex_trylock(m) ) ; ); 
    if ( result ) { result = TPI_ERROR_LOCK ; }
#else
    result = pthread_mutex_lock(m) ? TPI_ERROR_LOCK : 0 ;
#endif
  }

  return result ;
}

int TPI_Trylock( TPI_ThreadPool pool , int i )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result ) {
    result = 0 <= i && i < pool->m_mutex_size ? 0 : TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    result = pthread_mutex_trylock( pool->m_mutex + i );
    if ( EBUSY == result ) { result = TPI_LOCK_BUSY ; }
    else if ( result ) { result = TPI_ERROR_LOCK ; }
  }

  return result ;
}

int TPI_Unlock( TPI_ThreadPool pool , int i )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result ) {
    result = 0 <= i && i < pool->m_mutex_size ? 0 : TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    result = pthread_mutex_unlock( pool->m_mutex + i ) ? TPI_ERROR_LOCK : 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* Call the routine then change to NULL to indicate completed */

static void * local_thread_pool_driver( void * arg )
{
  ThreadData * const td = (ThreadData*)( arg );

  for ( int working = 1 ; working ; ) {

    if ( ! pthread_mutex_trylock( & td->m_lock ) ) {

      TPI_parallel_subprogram routine  = td->m_routine ;
      void                  * argument = td->m_argument ;

      /* When called by the driver, reset just the routine and argument */
      td->m_routine  = NULL ;
      td->m_argument = NULL ;

      if ( NULL != routine ) {
        (*routine)( argument, td->m_pool, td->m_size, td->m_rank );
      }

      if ( NULL == td->m_pool ) {
        working = 0 ;
      }

      pthread_mutex_unlock( & td->m_lock );
    }
  }

  return NULL ;
}

static void local_thread_pool_run_root( ThreadData * const td ,
                                        ThreadData * const td_root )
{
  TPI_ThreadPool          pool     = td->m_pool ;
  int                     size     = td->m_size ;
  int                     rank     = td->m_rank ;
  TPI_parallel_subprogram routine  = td->m_routine ;
  void                  * argument = td->m_argument ;

  td->m_routine  = td_root->m_routine  ;
  td->m_argument = td_root->m_argument ;
  td->m_pool     = td_root->m_pool ;
  td->m_size     = td_root->m_size ;
  td->m_rank     = td_root->m_rank ;

  if ( NULL != routine ) {
    (*routine)(argument,pool,size,rank);
  }
}

/*--------------------------------------------------------------------*/

int TPI_Run( TPI_ThreadPool pool ,
             TPI_parallel_subprogram routine ,
             void * routine_data )
{
  int result = 0 ;

  if ( NULL != routine && ! ( result = local_thread_pool_lock( pool ) ) ) {

    const int size = pool->m_end - pool->m_begin ;

    ThreadData tmp[ size ];
    int rank ;

    for ( rank = 0 ; rank < size ; ++rank ) {

      ThreadData * const td = pool->m_begin + rank ;
      ThreadData * const td_tmp = tmp + rank ;

      td_tmp->m_size     = td->m_size ;
      td_tmp->m_rank     = td->m_rank ;
      td_tmp->m_pool     = td->m_pool ;
      td_tmp->m_routine  = td->m_routine ;
      td_tmp->m_argument = td->m_argument ;

      td->m_size     = size ;
      td->m_rank     = rank ;
      td->m_pool     = pool ;
      td->m_routine  = routine ;
      td->m_argument = routine_data ;
    }

    for ( rank = 1 ; rank < size ; ++rank ) {
      pthread_mutex_unlock( & pool->m_begin[rank].m_lock );
    }

    /* Participate in the work */
    local_thread_pool_run_root( pool->m_begin , tmp );

    /* Wait for all threads to run and re-acquire locks */

    for ( int waiting = 1 ; waiting ; ) {

      waiting = 0 ;

      for ( rank = 1 ; rank < size ; ++rank ) {

        ThreadData * const td = pool->m_begin + rank ;
        ThreadData * const td_tmp = tmp + rank ;

        if ( td_tmp->m_pool != td->m_pool ) {
          if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
            local_thread_pool_run_root( td , td_tmp );
          }
        }
        if ( td_tmp->m_pool != td->m_pool ) { waiting = 1 ; }
      }
    }

    local_thread_pool_unlock( pool );
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Split( TPI_ThreadPool parent ,
               const int number ,
               const int sizes[] ,
               TPI_parallel_subprogram routine[] ,
               void * routine_data[] )
{
  int parent_size ;
  int result = TPI_Pool_size( parent , & parent_size );

  {
    int child_size = 0 ;

    for ( int i = 0 ; i < number && ! result ; ++i ) {
      if ( sizes[i] < 0 ) { result = TPI_ERROR_SIZE ; }
      child_size += sizes[i] ;
    }

    result = parent_size == child_size ? 0 : TPI_ERROR_SIZE ;
  }

  if ( ! result ) { result = local_thread_pool_lock( parent ); }

  if ( ! result ) {
    ThreadPool children[ number ];
    ThreadData tmp[ number ];
    int child , offset ;

    for ( offset = child = 0 ; child < number ; ++child ) {

      ThreadData * const td = parent->m_begin + offset ;

      if ( NULL != routine[child] && sizes[child] ) {

        ThreadData * const td_tmp = tmp + child ;

        td_tmp->m_size     = td->m_size ;
        td_tmp->m_rank     = td->m_rank ;
        td_tmp->m_pool     = td->m_pool ;
        td_tmp->m_routine  = td->m_routine ;
        td_tmp->m_argument = td->m_argument ;

        td->m_size     = number ;
        td->m_rank     = child ;
        td->m_pool     = children + child ;
        td->m_routine  = routine[child] ;
        td->m_argument = routine_data[child] ;

        children[child].m_begin = td ;
        children[child].m_end   = td + sizes[child] ;
        children[child].m_mutex_size = 0 ;

        result = local_thread_pool_lock_create( children + child );
      }

      offset += sizes[child] ;
    }

    for ( child = 1 ; child < number ; ++child ) {
      if ( NULL != routine[child] && sizes[child] ) {
        pthread_mutex_unlock( & children[child].m_begin->m_lock );
      }
    }

    if ( NULL != routine[0] ) {
      local_thread_pool_run_root( children->m_begin , tmp );
      TPI_Lock_allocation( children , -1 );
      local_thread_pool_lock_destroy( children );
    }

    for ( int waiting = 1 ; waiting ; ) {

      waiting = 0 ;

      for ( child = 1 ; child < number ; ++child ) {
        if ( NULL != routine[child] && sizes[child] ) { /* it was unlocked */

          ThreadData * const td     = children[child].m_begin ;
          ThreadData * const td_tmp = tmp + child ;

          if ( td_tmp->m_pool != td->m_pool ) { /* haven't processed it yet */
            if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
              local_thread_pool_run_root( td , td_tmp );
            }
            if ( td_tmp->m_pool != td->m_pool ) { waiting = 1 ; }
          }

          TPI_Lock_allocation( children + child , -1 );
          local_thread_pool_lock_destroy( children + child );
        }
      }
    }

    local_thread_pool_unlock( parent );
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void local_thread_pool_destroy_locks( int n , TPI_ThreadPool pool )
{
  while ( --n ) {
    ThreadData * const td = pool->m_begin + n ;
    pthread_mutex_unlock(  & td->m_lock );
    pthread_mutex_destroy( & td->m_lock );
  }
}

static int local_thread_pool_create_locks( int n , TPI_ThreadPool pool )
{
  int result = 0 ;
  int n_lock = 0 ;

  for ( n_lock = 1 ; n_lock < n && ! result ; ) {
    ThreadData * const td = pool->m_begin + n_lock ;
    if ( pthread_mutex_init( & td->m_lock , NULL ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else if ( pthread_mutex_lock( & td->m_lock ) ) {
      pthread_mutex_destroy( & td->m_lock );
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      ++n_lock ;
    }
  }

  if ( result ) {
    local_thread_pool_destroy_locks( n_lock , pool );
  }

  return result ;
}

static void local_thread_pool_verify_task(
  void * arg , TPI_ThreadPool pool , int size , int rank )
{
  ThreadData * const td = (ThreadData *) arg ;

  if ( NULL == td || NULL == pool ||
       size != pool->m_end - pool->m_begin ||
       rank != td - pool->m_begin ||
       pool != td->m_pool ||
       td->m_thread   != pthread_self() ) {
    td->m_rank = -1 ;
  }
}

static int local_thread_pool_verify_start( int n , TPI_ThreadPool pool )
{
  int result = 0 ;

  for ( int rank = 1 ; rank < n ; ++rank ) {
    ThreadData * td = pool->m_begin + rank ;

    td->m_pool = pool ;
    td->m_size = n ;
    td->m_rank = rank ;
    td->m_routine = & local_thread_pool_verify_task ;
    td->m_argument = td ;

    pthread_mutex_unlock( & td->m_lock );
  }

  for ( int waiting = 1 ; waiting ; ) {
    waiting = 0 ;

    for ( int rank = 1 ; rank < n ; ++rank ) {
      ThreadData * td = pool->m_begin + rank ;

      if ( NULL != td->m_pool ) {

        if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
          if ( NULL != td->m_routine ) { /* Has not run yet */
            pthread_mutex_unlock( & td->m_lock );
          }
          else { /* Has run */
            if ( -1 == td->m_rank ) { result = TPI_ERROR_INTERNAL ; }
            td->m_pool = NULL ;
            td->m_size = -1 ;
            td->m_rank = -1 ;
          }
        }

        if ( NULL != td->m_pool ) { waiting = 1 ; }
      }
    }
  }

  return result ;
}

static void local_thread_pool_terminate_threads( int n , TPI_ThreadPool pool )
{
  for ( int i = 1 ; i < n ; ++i ) {
    ThreadData * const td = pool->m_begin + i ;
    td->m_size     = 0 ;
    td->m_rank     = 0 ;
    td->m_pool     = NULL ;
    td->m_routine  = NULL ;
    td->m_argument = td ;
    pthread_mutex_unlock( & td->m_lock );
  }

  for ( int waiting = 1 ; waiting ; ) {
    waiting = 0 ;
    for ( int i = 1 ; i < n ; ++i ) {
      ThreadData * const td = pool->m_begin + i ;

      if ( ! td->m_size ) {
        if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
          if ( NULL != td->m_argument ) {
            pthread_mutex_unlock( & td->m_lock );
          }
          else {
            td->m_size = -1 ;
            td->m_rank = -1 ;
          }
        }
        if ( ! td->m_size ) { waiting = 1 ; }
      }
    }
  }
}

static int local_thread_pool_create_threads( int n , TPI_ThreadPool pool )
{
  int result = 0 ;
  int n_thread ;
  pthread_attr_t thread_attr ;

  for ( n_thread = 0 ; n_thread < n ; ++n_thread ) {

    ThreadData * const td = pool->m_begin + n_thread ;

    td->m_size     = -1 ;
    td->m_rank     = -1 ;
    td->m_pool     = NULL ;
    td->m_routine  = NULL ;
    td->m_argument = NULL ;
  }

  /* pthread_setconcurrency( n ); */

  if ( pthread_attr_init( & thread_attr ) ) {
    result = TPI_ERROR_INTERNAL ;
  }

  if ( ! result ) {
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

/* GNU runtime rejects this attribute:
    pthread_attr_setschedpolicy( & thread_attr, SCHED_FIFO );
*/

  }

  for ( n_thread = 1 ; n_thread < n && ! result ; ) {

    ThreadData * const td = pool->m_begin + n_thread ;

    const int p_result =
      pthread_create( & td->m_thread ,
                      & thread_attr ,
                      & local_thread_pool_driver ,
                      (void *) td );

    if ( p_result ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      ++n_thread ;
    }
  }

  pthread_attr_destroy( & thread_attr );

  if ( ! result ) {
    pool->m_begin->m_thread = pthread_self();
    pool->m_end = pool->m_begin + n_thread ;
    result = local_thread_pool_verify_start( n_thread , pool );
  }

  if ( result ) {
    pool->m_end = pool->m_begin ;
    local_thread_pool_terminate_threads( n_thread , pool );
  }

  return result ;
}

int TPI_Init( int n , TPI_ThreadPool * root )
{
  static ThreadData threads[ MAX_THREADS ];
  static int ultimate_init = 1 ;

  int result = 0 < n && n < MAX_THREADS ? 0 : TPI_ERROR_SIZE ;

  if ( ! result && NULL == root ) { result = TPI_ERROR_NULL ; }

  if ( ! result ) {

    TPI_ThreadPool pool = local_thread_pool_root() ;

    if ( ultimate_init ) {
      pool->m_begin = threads ;
      pool->m_end   = threads ;
      pool->m_mutex_size = 0 ;
      ultimate_init = 0 ;
    }

    if ( pool->m_begin != pool->m_end ) {
      result = TPI_ERROR_ACTIVE ;
    }

    if ( ! result ) { result = local_thread_pool_lock_create( pool ); }

    if ( ! result ) {
      result = local_thread_pool_create_locks( n , pool );
    }

    if ( ! result ) {
      /* Have locks */

      result = local_thread_pool_create_threads( n , pool );

      if ( result ) {
        local_thread_pool_destroy_locks( n , pool );
      }
    }
    *root = result ? NULL : pool ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  TPI_ThreadPool pool = local_thread_pool_root() ;

  int result = local_thread_pool_lock( pool );

  if ( ! result ) {

    local_thread_pool_unlock( pool );
    local_thread_pool_lock_destroy( pool );

    const int n = pool->m_end - pool->m_begin ;

    local_thread_pool_terminate_threads( n , pool );
    local_thread_pool_destroy_locks( n , pool );

    TPI_Lock_allocation( pool , -1 );

    pool->m_end = pool->m_begin ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef struct TPI_ThreadPool_Private {
  TPI_ThreadPool m_root ;
  int            m_size ;
  int            m_lock_size ;
  int            m_lock[ MAX_LOCK + 1 ];
} ThreadPool ;

/*--------------------------------------------------------------------*/

static TPI_ThreadPool local_thread_pool_root()
{
  static ThreadPool root_pool = { NULL , 0 , 0 };
  return & root_pool ;
}

int TPI_Pool_size( TPI_ThreadPool pool , int * size )
{
  int result = NULL != pool ? 0 : TPI_ERROR_NULL ;

  if ( ! result && pool->m_root != local_thread_pool_root() ) {
    result = TPI_ERROR_POOL ;
  }

  if ( ! result && size ) { *size = pool->m_size ; }

  return result ;
}

int TPI_Lock_size( TPI_ThreadPool pool , int * size )
{
  int result = TPI_Pool_size( pool , NULL );
  if ( ! result && size ) { *size = pool->m_lock_size ; }
  return result ;
}

int TPI_Lock_allocation( TPI_ThreadPool pool , int number )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result && pool->m_lock[ MAX_LOCK ] ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result ) {
    if ( MAX_LOCK < number ) {
      result = TPI_ERROR_SIZE ;
    }
    else if ( number < 0 ) {
      for ( int i = 0 ; i < pool->m_lock_size ; ++i ) {
        pool->m_lock[i] = 0 ;
      }
      pool->m_lock_size = 0 ;
    }
    else if ( pool->m_lock_size < number ) {
      for ( int i = pool->m_lock_size ; i < number ; ++i ) {
        pool->m_lock[i] = 0 ;
      }
      pool->m_lock_size = number ;
    }
  }
  return result ;
}

int TPI_Run( TPI_ThreadPool pool ,
             TPI_parallel_subprogram routine ,
             void * routine_data )
{
  int size ;
  int result = TPI_Pool_size( pool , & size );

  if ( ! result && pool->m_lock[ MAX_LOCK ] ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result && NULL != routine ) {
    pool->m_lock[ MAX_LOCK ] = 1 ;
    for ( int rank = 0 ; rank < size ; ++rank ) {
      (*routine)( routine_data , pool , size , rank );
    }
    pool->m_lock[ MAX_LOCK ] = 0 ;
  }
  return result ;
}

int TPI_Lock( TPI_ThreadPool pool , int i )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result && pool->m_lock_size <= i ) { result = TPI_ERROR_SIZE ; }

  if ( ! result && pool->m_lock[i] ) { result = TPI_ERROR_LOCK ; }

  if ( ! result ) { pool->m_lock[i] = 1 ; }

  return result ;
}

int TPI_Split( TPI_ThreadPool parent ,
               const int number ,
               const int sizes[] ,
               TPI_parallel_subprogram routine[] ,
               void * routine_data[] )
{
  int parent_size ;
  int result = TPI_Pool_size( parent , & parent_size );

  {
    int child_size = 0 ;

    for ( int i = 0 ; i < number && ! result ; ++i ) {
      if ( sizes[i] < 0 ) { result = TPI_ERROR_SIZE ; }
      child_size += sizes[i] ;
    }

    result = parent_size == child_size ? 0 : TPI_ERROR_SIZE ;
  }

  if ( ! result && parent->m_lock[MAX_LOCK] ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result ) {
    TPI_ThreadPool root = local_thread_pool_root();
    ThreadPool children[ number ];

    for ( int child = 0 ; child < number ; ++child ) {
      children[child].m_root = root ;
      children[child].m_size = sizes[child] ;
      children[child].m_lock_size = 0 ;
      children[child].m_lock[ MAX_LOCK ] = 0 ;

      if ( NULL != routine[child] && sizes[child] ) {
        (*routine[child])( routine_data[child] , children + child ,
                           number , child );
      }
    }
  }

  return result ;
}

int TPI_Trylock( TPI_ThreadPool pool , int i )
{ return TPI_Lock( pool , i ); }

int TPI_Unlock( TPI_ThreadPool pool , int i )
{
  int result = TPI_Pool_size( pool , NULL );

  if ( ! result && pool->m_lock_size <= i ) { result = TPI_ERROR_SIZE ; }

  if ( ! result && ! pool->m_lock[i] ) { result = TPI_ERROR_LOCK ; }

  if ( ! result ) { pool->m_lock[i] = 0 ; }

  return result ;
}

int TPI_Init( int n , TPI_ThreadPool * root )
{
  TPI_ThreadPool pool = local_thread_pool_root() ;

  int result = 0 < n && n < MAX_THREADS ? 0 : TPI_ERROR_SIZE ;

  if ( ! result ) { result = root != NULL ? 0 : TPI_ERROR_NULL ; }

  if ( ! result && pool->m_size ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result ) {
    pool->m_root = pool ;
    pool->m_size = n ;
    pool->m_lock_size = 0 ;
    pool->m_lock[ MAX_LOCK ] = 0 ;
    *root = pool ;
  }

  return result ;
}

int TPI_Finalize()
{
  TPI_ThreadPool pool = local_thread_pool_root();

  int result = pool->m_lock[ MAX_LOCK ] ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {
    TPI_Lock_allocation( pool , -1 );
    pool->m_root = NULL ;
    pool->m_size = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#endif

