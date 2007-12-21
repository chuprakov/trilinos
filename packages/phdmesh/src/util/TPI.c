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

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#ifdef NO_PTHREADS_AVAILABLE

int TPI_Init( int , TPI_ThreadPool * root )
{ if ( root ) { *root = NULL ; } return 0 ; }

int TPI_Finalize() {}

int TPI_Pool_size( TPI_ThreadPool , int * size )
{ if ( size ) { *size = 1 ; } return 0 ; }

int TPI_Lock_size( TPI_ThreadPool , int * size )
{ if ( size ) { *size = 0 ; } return 0 ; }

int TPI_Run( TPI_ThreadPool ,
             TPI_parallel_subprogram prog , void * arg , int )
{ if ( prog ) { (*prog)( arg , NULL , 0 ); } return 0 ; }

int TPI_Lock( TPI_ThreadPool , int )
{ return 0 ; }

int TPI_Unlock( TPI_ThreadPool , int )
{ return 0 ; }

int TPI_Trylock( TPI_ThreadPool , int )
{ return 0 ; }

int TPI_Split( TPI_ThreadPool ,
               const int n ,
               const int sizes[] ,
               TPI_parallel_subprogram prog[] ,
               void * arg[] ,
               int )
{
  int result = 1 == n && sizes[0] == 1 ? 0 : TPI_ERROR_SIZE ;

  if ( ! result && prog[0] ) { (*(prog[0]))( arg[0] , NULL , 0 ); }

  return result ;
}

int TPI_Split_lock_size( TPI_ThreadPool , int * size )
{ if ( size ) { *size = 0 ; } return 0 ; }

int TPI_Split_lock(      TPI_ThreadPool , int ) { return 0 ; }
int TPI_Split_trylock(   TPI_ThreadPool , int ) { return 0 ; }
int TPI_Split_unlock(    TPI_ThreadPool , int ) { return 0 ; }

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else

#include <unistd.h>
#include <errno.h>
#include <pthread.h>

typedef struct ThreadDataStruct {
  pthread_t               m_thread ;
  pthread_mutex_t         m_lock ;
  TPI_ThreadPool          m_pool ;
  int                     m_rank ;
  TPI_parallel_subprogram m_routine ;
  void                  * m_argument ;
} ThreadData ;

typedef struct TPI_ThreadPool_Private {
  pthread_mutex_t   m_active ;
  pthread_mutex_t * m_mutex ;
  int               m_mutex_size ;
  ThreadData      * m_begin ;
  ThreadData      * m_end ;
  TPI_ThreadPool    m_parent ;
  void           (* m_check )();
} ThreadPool ;

static void local_thread_pool_check() {}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Pool_size( TPI_ThreadPool pool , int * size )
{
  int result = 0 ;

  if ( NULL != pool ) {
    if ( & local_thread_pool_check != pool->m_check ) {
      result = TPI_ERROR_POOL ;
    }
    else if ( size ) {
      *size = pool->m_end - pool->m_begin ;
    }
  }
  else if ( size ) {
    *size = 1 ;
  }

  return result ;
}

int TPI_Lock_size( TPI_ThreadPool pool , int * size )
{
  int result = 0 ;

  if ( NULL != pool ) {
    if ( & local_thread_pool_check != pool->m_check ) {
      result = TPI_ERROR_POOL ;
    }
    else if ( size ) {
      *size = pool->m_mutex_size ;
    }
  }
  else if ( size ) {
    *size = 0 ;
  }

  return result ;
}

int TPI_Split_lock_size( TPI_ThreadPool pool , int * size )
{
  int result = 0 ;

  if ( NULL != pool ) {
    if ( & local_thread_pool_check != pool->m_check ) {
      result = TPI_ERROR_POOL ;
    }
    else if ( size ) {
      *size = pool->m_parent ? pool->m_parent->m_mutex_size : 0 ;
    }
  }
  else if ( size ) {
    *size = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static
int local_thread_pool_mutex( TPI_ThreadPool pool , int i ,
                             pthread_mutex_t ** m )
{
  int result = NULL != pool ? 0 : TPI_ERROR_NULL ;

  if ( ! result && & local_thread_pool_check != pool->m_check ) {
    result = TPI_ERROR_POOL ;
  }

  if ( ! result && ( i < 0 || pool->m_mutex_size <= i ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) { *m = pool->m_mutex + i ; }

  return result ;
}

static
int local_thread_pool_split_mutex( TPI_ThreadPool pool , int i ,
                                   pthread_mutex_t ** m )
{
  int result = NULL != pool ? 0 : TPI_ERROR_NULL ;

  if ( ! result && & local_thread_pool_check != pool->m_check ) {
    result = TPI_ERROR_POOL ;
  }

  if ( ! result && ( NULL == pool->m_parent || i < 0 ||
                     pool->m_parent->m_mutex_size <= i ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) { *m = pool->m_parent->m_mutex + i ; }

  return result ;
}

static
int local_thread_pool_lock( pthread_mutex_t * m )
{
  int result ;
#if 1
  while ( EBUSY == ( result = pthread_mutex_trylock(m) ) ); 
  if ( result ) { result = TPI_ERROR_LOCK ; }
#else
  result = pthread_mutex_lock(m) ? TPI_ERROR_LOCK : 0 ;
#endif
  return result ;
}

static
int local_thread_pool_trylock( pthread_mutex_t * m )
{
  int result = pthread_mutex_trylock( m );
  if ( EBUSY == result ) { result = TPI_LOCK_BUSY ; }
  else if ( result ) { result = TPI_ERROR_LOCK ; }
  return result ;
}

int TPI_Lock( TPI_ThreadPool pool , int i )
{
  pthread_mutex_t * m ;

  int result = local_thread_pool_mutex( pool , i , & m );

  if ( ! result ) { result = local_thread_pool_lock( m ); }

  return result ;
}

int TPI_Split_lock( TPI_ThreadPool pool , int i )
{
  pthread_mutex_t * m ;

  int result = local_thread_pool_split_mutex( pool , i , & m );

  if ( ! result ) { result = local_thread_pool_lock( m ); }

  return result ;
}

int TPI_Trylock( TPI_ThreadPool pool , int i )
{
  pthread_mutex_t * m ;

  int result = local_thread_pool_mutex( pool , i , & m );

  if ( ! result ) { result = local_thread_pool_trylock( m ); }

  return result ;
}

int TPI_Split_trylock( TPI_ThreadPool pool , int i )
{
  pthread_mutex_t * m ;

  int result = local_thread_pool_split_mutex( pool , i , & m );

  if ( ! result ) { result = local_thread_pool_trylock( m ); }

  return result ;
}

int TPI_Unlock( TPI_ThreadPool pool , int i )
{
  pthread_mutex_t * m ;

  int result = local_thread_pool_mutex( pool , i , & m );

  if ( ! result ) { result = pthread_mutex_unlock( m ) ? TPI_ERROR_LOCK : 0 ; }

  return result ;
}

int TPI_Split_unlock( TPI_ThreadPool pool , int i )
{
  pthread_mutex_t * m ;

  int result = local_thread_pool_split_mutex( pool , i , & m );

  if ( ! result ) { result = pthread_mutex_unlock( m ) ? TPI_ERROR_LOCK : 0 ; }

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
        (*routine)( argument, td->m_pool, td->m_rank );
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
  int                     rank     = td->m_rank ;
  TPI_parallel_subprogram routine  = td->m_routine ;
  void                  * argument = td->m_argument ;

  td->m_routine  = td_root->m_routine  ;
  td->m_argument = td_root->m_argument ;
  td->m_pool     = td_root->m_pool ;
  td->m_rank     = td_root->m_rank ;

  if ( NULL != routine ) {
    (*routine)(argument,pool,rank);
  }
}

/*--------------------------------------------------------------------*/

static void local_thread_pool_run( TPI_ThreadPool pool ,
                                   TPI_parallel_subprogram routine ,
                                   void * routine_data )
{
  const int size = pool->m_end - pool->m_begin ;

  ThreadData tmp[ size ];

  for ( int rank = 0 ; rank < size ; ++rank ) {

    ThreadData * const td = pool->m_begin + rank ;
    ThreadData * const td_tmp = tmp + rank ;

    td_tmp->m_rank     = td->m_rank ;
    td_tmp->m_pool     = td->m_pool ;
    td_tmp->m_routine  = td->m_routine ;
    td_tmp->m_argument = td->m_argument ;

    td->m_rank     = rank ;
    td->m_pool     = pool ;
    td->m_routine  = routine ;
    td->m_argument = routine_data ;
  }

  for ( int rank = 1 ; rank < size ; ++rank ) {
    pthread_mutex_unlock( & pool->m_begin[rank].m_lock );
  }

  /* Participate in the work */
  local_thread_pool_run_root( pool->m_begin , tmp );

  /* Wait for all threads to run and re-acquire locks */

  for ( int waiting = 1 ; waiting ; ) {

    waiting = 0 ;

    for ( int rank = 1 ; rank < size ; ++rank ) {

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
}

int TPI_Run( TPI_ThreadPool pool ,
             TPI_parallel_subprogram routine ,
             void * routine_data ,
             int number_locks )
{
  int result = number_locks < 0 ? TPI_ERROR_SIZE : 0 ;

  if ( ! result && NULL != routine ) {

    pthread_mutex_t locks[ number_locks ];

    int nlocks = 0 ;

    while ( nlocks < number_locks && ! result ) {
      if ( ! pthread_mutex_init( locks + nlocks , NULL ) ) {
        ++nlocks ;
      }
      else {
        result = TPI_ERROR_INTERNAL ;
      }
    }

    if ( ! result ) {

      if ( NULL != pool ) {

        if ( & local_thread_pool_check != pool->m_check ) {
          result = TPI_ERROR_POOL ;
        }

        if ( ! result && pthread_mutex_trylock( & pool->m_active ) ) {
          result = TPI_ERROR_ACTIVE ;
        }

        if ( ! result ) {
 
          pool->m_mutex = locks ;
          pool->m_mutex_size = number_locks ;

          local_thread_pool_run( pool , routine , routine_data );

          pool->m_mutex = NULL ;
          pool->m_mutex_size = 0 ;

          pthread_mutex_unlock( & pool->m_active );
        }
      }
      else {
        ThreadData dtmp = { pthread_self() , PTHREAD_MUTEX_INITIALIZER ,
                            NULL , -1 , NULL , NULL };

        ThreadPool ptmp = { PTHREAD_MUTEX_INITIALIZER,
                            locks , number_locks ,
                            & dtmp, ( & dtmp ) + 1 ,
                            NULL , & local_thread_pool_check };

        dtmp.m_pool = & ptmp ;

        if ( ! pthread_mutex_trylock( & ptmp.m_active ) ) {
 
          local_thread_pool_run( & ptmp , routine , routine_data );

          pthread_mutex_unlock(  & ptmp.m_active );
          pthread_mutex_destroy( & ptmp.m_active );
        }
        else {
          result = TPI_ERROR_INTERNAL ;
        }
      }
    }

    while ( nlocks-- ) {
      pthread_mutex_destroy( locks + nlocks );
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int local_thread_pool_split_run( TPI_ThreadPool parent ,
                                 const int number ,
                                 const int sizes[] ,
                                 TPI_parallel_subprogram routine[] ,
                                 void * routine_data[] )
{
  ThreadPool children[ number ];
  ThreadData tmp[ number ];
  int child , offset ;
  int result = 0 ;

  for ( offset = child = 0 ; child < number && ! result ; ++child ) {

    ThreadData * const td = parent->m_begin + offset ;

    if ( NULL != routine[child] && sizes[child] ) {

      ThreadData * const td_tmp = tmp + child ;

      td_tmp->m_rank     = td->m_rank ;
      td_tmp->m_pool     = td->m_pool ;
      td_tmp->m_routine  = td->m_routine ;
      td_tmp->m_argument = td->m_argument ;

      td->m_rank     = 0 ;
      td->m_pool     = children + child ;
      td->m_routine  = routine[child] ;
      td->m_argument = routine_data[child] ;

      children[child].m_mutex  = NULL ;
      children[child].m_mutex_size = 0 ;
      children[child].m_begin  = td ;
      children[child].m_end    = td + sizes[child] ;
      children[child].m_parent = parent ;
      children[child].m_check  = & local_thread_pool_check ;

      if ( pthread_mutex_init( & children[child].m_active , NULL ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
    }

    offset += sizes[child] ;
  }

  if ( result ) {
    while ( child-- ) {
      pthread_mutex_destroy( & children[child].m_active );
    }
  }
  else {

    for ( child = 1 ; child < number ; ++child ) {
      if ( NULL != routine[child] && sizes[child] ) {
        pthread_mutex_unlock( & children[child].m_begin->m_lock );
      }
    }

    if ( NULL != routine[0] ) {
      local_thread_pool_run_root( children->m_begin , tmp );
      pthread_mutex_destroy( & children->m_active );
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

          pthread_mutex_destroy( & children[child].m_active );
        }
      }
    }
  }

  return result ;
}

int TPI_Split( TPI_ThreadPool parent ,
               const int number ,
               const int sizes[] ,
               TPI_parallel_subprogram routine[] ,
               void * routine_data[] ,
               int number_locks )
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

  if ( ! result && pthread_mutex_trylock( & parent->m_active ) ) {
    result = TPI_ERROR_ACTIVE ;
  }

  if ( ! result ) {

    pthread_mutex_t locks[ number_locks ];

    for ( int i = 0 ; i < number_locks && ! result ; ++i ) {
      if ( pthread_mutex_init( locks + i , NULL ) ) {
        result = TPI_ERROR_INTERNAL ;
        while ( i-- ) { pthread_mutex_destroy( locks + i ); }
      }
    }

    if ( ! result ) {

      parent->m_mutex = locks ;
      parent->m_mutex_size = number_locks ;

      result = local_thread_pool_split_run( parent , number , sizes ,
                                            routine , routine_data );

      parent->m_mutex = NULL ;
      parent->m_mutex_size = 0 ;

      for ( int i = 0 ; i < number_locks ; ++i ) {
        pthread_mutex_destroy( locks + i );
      }
    }

    pthread_mutex_unlock( & parent->m_active );
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
  void * arg , TPI_ThreadPool pool , int rank )
{
  ThreadData * const td = (ThreadData *) arg ;

  if ( NULL == td ||
       NULL == pool ||
       pool         != td->m_pool ||
       rank         != td - pool->m_begin ||
       td->m_thread != pthread_self() ) {
    td->m_rank = -1 ;
  }
}

static int local_thread_pool_verify_start( int n , TPI_ThreadPool pool )
{
  int result = 0 ;

  for ( int rank = 1 ; rank < n ; ++rank ) {
    ThreadData * td = pool->m_begin + rank ;

    td->m_pool = pool ;
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

      if ( ! td->m_rank ) {
        if ( ! pthread_mutex_trylock( & td->m_lock ) ) {
          if ( NULL != td->m_argument ) {
            pthread_mutex_unlock( & td->m_lock );
          }
          else {
            td->m_rank = -1 ;
          }
        }
        if ( ! td->m_rank ) { waiting = 1 ; }
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

static int local_thread_pool_root( int n , TPI_ThreadPool * pool )
{
  enum { MAX_THREADS = 1024 };

  static ThreadData threads[ MAX_THREADS ];

  static ThreadPool root = {
      PTHREAD_MUTEX_INITIALIZER ,
      NULL , 0 ,
      threads , threads ,
      NULL ,
      & local_thread_pool_check };

  int result = threads == root.m_begin ? 0 : TPI_ERROR_INTERNAL ;

  if ( ! result && pthread_mutex_trylock( & root.m_active ) ) {
    result = TPI_ERROR_ACTIVE ;
  }

  if ( ! result ) {

    if ( pool ) { /* Initialize */

      result = 0 < n && n < MAX_THREADS ? 0 : TPI_ERROR_SIZE ;

      if ( ! result && threads != root.m_end ) { result = TPI_ERROR_ACTIVE ; }

      if ( ! result ) {
        result = local_thread_pool_create_locks( n , & root );
      }

      if ( ! result ) { /* Have locks */

        result = local_thread_pool_create_threads( n , & root );

        if ( result ) {
          local_thread_pool_destroy_locks( n , & root );
        }
      }
      *pool = result ? NULL : & root ;
    }
    else if ( root.m_begin < root.m_end ) { /* Shut down */

      n = root.m_end - root.m_begin ;

      local_thread_pool_terminate_threads( n , & root );
      local_thread_pool_destroy_locks(     n , & root );

      root.m_end = root.m_begin ;
    }

    pthread_mutex_unlock( & root.m_active );
  }

  return result ;
}

int TPI_Init( int n , TPI_ThreadPool * pool )
{
  int result = NULL != pool ? 0 : TPI_ERROR_NULL ;

  if ( ! result ) { result = local_thread_pool_root( n , pool ); }

  return result ;
}

int TPI_Finalize()
{
  return local_thread_pool_root( 0 , NULL );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#endif

