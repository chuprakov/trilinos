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
#include <util/taskpool.h>

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

enum LockControl { LOCK , UNLOCK , RESIZE , CLEAR };

static
void phdmesh_taskpool_lock_control(
  enum LockControl control ,
  unsigned n ,
  const char * const error_info );

static
int phdmesh_taskpool_control(
  unsigned                 number ,
  phdmesh_taskpool_routine routine ,
  void                   * routine_data );

void phdmesh_taskpool_lock( unsigned n , const char * const error_info )
{ phdmesh_taskpool_lock_control( LOCK , n , error_info ); }

void phdmesh_taskpool_unlock( unsigned n , const char * const error_info )
{ phdmesh_taskpool_lock_control( UNLOCK , n , error_info ); }

int phdmesh_taskpool_resize( unsigned n )
{ return phdmesh_taskpool_control( n , NULL , NULL ); }

int phdmesh_taskpool_run(
  phdmesh_taskpool_routine routine ,
  void * routine_data ,
  unsigned      number )
{
  return NULL != routine ?
         phdmesh_taskpool_control( number, routine, routine_data) : 0 ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if defined( PHDMESH_HAS_PTHREADS )

#include <pthread.h>

static void phdmesh_private_lock( pthread_mutex_t * m )
{
#if 1
   while ( pthread_mutex_trylock( m ) );
#else
   pthread_mutex_lock( m );
#endif
}

static
void phdmesh_taskpool_lock_control(
  enum LockControl control , unsigned n , const char * const error_info )
{
  enum { MAX_LOCKS = 512 };

  static pthread_mutex_t mutex[ MAX_LOCKS ];
  static unsigned number = 0 ;

  int ok = n < MAX_LOCKS ;

  if ( ok ) {
    switch( control ) {
    case LOCK :
      ok = n < number ;
      if ( ok ) { phdmesh_private_lock( mutex + n ); }
      break ;
    case UNLOCK :
      ok = n < number ;
      if ( ok ) { pthread_mutex_unlock( mutex + n ); }
      break ;
    case RESIZE :
      while ( number < n ) {
        if ( pthread_mutex_init( mutex + number , NULL ) ) { ok = 0 ; }
        ++number ;
      }
      break ;
    case CLEAR :
      while ( number ) {
        --number ;
        if ( pthread_mutex_destroy( mutex + number ) ) { ok = 0 ; }
      }
      break ;
    }
  }

  if ( ! ok ) {
    static const char * const op_text[] = { "lock" , "unlock" , "run" };
    const char * const info = error_info ? error_info : "" ;

    fprintf( stderr , "phdmesh_taskpool_%s( %u , %s )[ %u of %u ] FAILED\n" ,
             op_text[ control ] , n , info , number , (unsigned) MAX_LOCKS );

    abort();
  }
}


typedef struct ThreadDataStruct {
  pthread_t                m_thread ;
  pthread_mutex_t          m_active ;
  pthread_mutex_t          m_sync ;
  phdmesh_taskpool_routine m_routine ;
  void                   * m_argument ;
  unsigned                 m_size ;
  unsigned                 m_rank ;
  int                      m_result ;
} ThreadData ;

static
void * phdmesh_taskpool_driver( void * arg )
{
  while ( NULL != arg ) {
    ThreadData * const t = (ThreadData*)( arg );

    phdmesh_private_lock( & t->m_active );

    /* Own non-mutex data when activate */

    if ( NULL == t->m_routine ) {
      arg = NULL ;
      t->m_result = 0 ;
    }
    else {
      t->m_result = (*t->m_routine)( t->m_argument , t->m_size , t->m_rank );
    }
    t->m_routine = NULL ;
    t->m_argument = NULL ;

    pthread_mutex_unlock( & t->m_active ); /* Done with activation */

    phdmesh_private_lock( & t->m_sync ); /* Waiting to sync */

    pthread_mutex_unlock( & t->m_sync );   /* Done sync, am now inactive */
  }

  return NULL ;
}

static
int phdmesh_taskpool_control(
  unsigned                 number ,
  phdmesh_taskpool_routine routine ,
  void                   * routine_data )
{
  enum { MAX_PTHREADS = 512 };
  static ThreadData threads[ MAX_PTHREADS ];
  static unsigned num_threads = 0 ;
  static int active = 0 ;

  int result = active ;

  if ( ! active ) {

    active = -1 ;

    /*----------------------------------------------------------------*/
    if ( NULL != routine ) {

      ThreadData * tptr[ num_threads ];

      phdmesh_taskpool_lock_control( RESIZE , number , NULL );

      for ( unsigned i = 0 ; i < num_threads ; ++i ) {
        ThreadData * const t = tptr[i] = threads + i ;
        t->m_routine  = routine ;
        t->m_argument = routine_data ;

        /* The completion rendezvous lock is always available */
        pthread_mutex_lock(   & t->m_sync );
        pthread_mutex_unlock( & t->m_active ); /* Unblock worker */
      }

      /* Participate in the work as p_rank == 0. */
      result = (*routine)( routine_data , num_threads + 1 , 0 );

      unsigned n = num_threads ;

      while ( n ) {
        for ( unsigned i = 0 ; i < num_threads ; ++i ) {
          ThreadData * const t = tptr[i] ;

          if ( NULL != t && ! pthread_mutex_trylock( & t->m_active ) ) {
            pthread_mutex_unlock( & t->m_sync );  /* Allow worker to loop */

            if ( NULL != t->m_routine ) {
              /* Worker thread did not run */
              t->m_result = (*t->m_routine)(t->m_argument,t->m_size,t->m_rank);
              t->m_routine = NULL ;
              t->m_argument = NULL ;
            }

            result |= t->m_result ;

            --n ;
          }
        }
      }
    }
    /*----------------------------------------------------------------*/
    else {

      if ( number ) {
        --number ; /* The 'main' is one of the threads */
      }
      else {
        phdmesh_taskpool_lock_control( CLEAR , 0 , NULL );
      }

      const unsigned n = number < MAX_PTHREADS ? number : MAX_PTHREADS ;

      while ( n < num_threads ) {
        --num_threads ;
        pthread_mutex_unlock( & threads[ num_threads ].m_active );
        pthread_mutex_destroy( & threads[ num_threads ].m_active );
        pthread_mutex_destroy( & threads[ num_threads ].m_sync );
      }

      if ( num_threads < n ) {

        /* pthread_setconcurrency( n + 1 ); */

        pthread_attr_t thread_attr ;
        pthread_attr_init( & thread_attr );
        pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
        pthread_attr_setschedpolicy( & thread_attr, SCHED_FIFO );
        pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED);

        while ( ! result && num_threads < n ) {

          ThreadData * const data = threads + num_threads ;

          data->m_routine  = NULL ;
          data->m_argument = NULL ;
          data->m_size     = 0 ;
          data->m_rank     = num_threads + 1 ;
          data->m_result   = 0 ;

          pthread_mutex_init( & data->m_sync ,   NULL );
          pthread_mutex_init( & data->m_active , NULL );

          pthread_mutex_lock( & data->m_active );

          result == pthread_create( & data->m_thread , & thread_attr ,
                                    & phdmesh_taskpool_driver ,
                                    (void *) data );
          if ( ! result ) { ++num_threads ; }
        }

        pthread_attr_destroy( & thread_attr );
      }

      for ( unsigned i = 0 ; i < n ; ++i ) {
        threads[i].m_size = n + 1 ;
      }
    }

    active = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/

#else /* ! defined( PHDMESH_HAS_PTHREADS ) */

static
void phdmesh_taskpool_lock_control(
  enum LockControl control , unsigned n , const char * const error_info )
{
  enum { MAX_LOCKS = 512 };

  static unsigned mutex[ MAX_LOCKS ];
  static unsigned number = 0 ;

  int ok = n < MAX_LOCKS ;

  if ( ok ) {
    switch( control ) {
    case LOCK :
      ok = n < number && ! mutex[n] ;
      if ( ok ) { mutex[n] = 1 ; }
      break ;
    case UNLOCK :
      ok = n < number && mutex[n] ;
      if ( ok ) { mutex[n] = 0 ; }
      break ;
    case RESIZE :
      while ( number < n ) { mutex[ number ] = 0 ; ++number ; }
      break ;
    case CLEAR :
      while ( number ) { --number ; mutex[ number ] = 0 ; }
      break ;
    }
  }

  if ( ! ok ) {
    static const char * const op_text[] = { "lock" , "unlock" , "run" };

    fprintf( stderr , "phdmesh_taskpool_%s( %u , %s )[ %u of %u ] FAILED\n" ,
             op_text[ control ] , n , error_info ,
             number , (unsigned) MAX_LOCKS );

    abort();
  }
}

static
int phdmesh_taskpool_control(
  unsigned      number ,
  phdmesh_taskpool_routine routine ,
  void * routine_data )
{
  static int active = 0 ;

  int result = active ;

  if ( ! active ) {
    active = -1 ;

    if ( NULL != routine ) {

      phdmesh_taskpool_lock_control( RESIZE , number , NULL );

      result = (*routine)( routine_data , 1 , 0 );
    }
    else {
      phdmesh_taskpool_lock_control( CLEAR , 0 , NULL );
    }

    active = 0 ;
  }

  return result ;
}

#endif

