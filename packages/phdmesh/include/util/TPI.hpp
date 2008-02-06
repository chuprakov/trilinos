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

#ifndef util_ThreadPool_hpp
#define util_ThreadPool_hpp

#include <string>
#include <stdexcept>

#include <util/TPI.h>

namespace TPI {

typedef TPI_ThreadPool ThreadPool ;

inline
int Pool_rank( ThreadPool pool , int & rank , int & size )
  { return TPI_Pool_rank( pool , & rank , & size ); }

inline
int Set_buffer_size( ThreadPool pool , int n )
  { return TPI_Set_buffer_size( pool , n ); }

inline
int Buffer( ThreadPool pool , void * & buf , int & size )
  { return TPI_Buffer( pool , & buf , & size ); }

inline
int Lock_size( ThreadPool pool , int & n )
  { return TPI_Lock_size( pool , & n ); }

inline
int Set_lock_size( ThreadPool pool , int n )
  { return TPI_Set_lock_size( pool , n ); }

inline
int Lock( ThreadPool pool , int n ) { return TPI_Lock( pool , n ); }

inline
int Unlock( ThreadPool pool , int n ) { return TPI_Unlock( pool , n ); }

inline
void Run( ThreadPool, void (*)(void*,ThreadPool), void * );

/** Run 'void Worker::method( TPI::ThreadPool )' */
template<class Worker>
inline
void Run( ThreadPool , Worker & , void (Worker::*)(ThreadPool) );

/** Lock guard to insure that a lock is released
 *  when control exists a block.
 *    {
 *      TPI::LockGuard local_lock( i );
 *    }
 */
template<class Worker>
class LockGuard {
private:
  LockGuard();
  LockGuard( const LockGuard & );
  LockGuard & operator = ( const LockGuard & );
  const ThreadPool m_pool ;
  const int        m_value ;
  const int        m_result ;
public:
  operator int() const { return m_result ; }

  explicit LockGuard( ThreadPool pool , unsigned i_lock )
    : m_pool( pool ), m_value( i_lock ), m_result( TPI_Lock(pool,i_lock) ) {}

  ~LockGuard() { TPI_Unlock( m_pool , m_value ); }
};

//----------------------------------------------------------------------

namespace {

template<class Worker>
class WorkerMethod {
private:
  WorkerMethod();
  WorkerMethod( const WorkerMethod & );
  WorkerMethod & operator = ( const WorkerMethod & );

public:

  typedef void (Worker::*Method)( ThreadPool );

  Worker & worker ;
  Method   method ;

  WorkerMethod( Worker & w , Method m ) : worker(w), method(m) {}

  void run_p( ThreadPool pool ) const
    { try { (worker.*method)(pool); } catch(...){} }

  static void run( void * arg , ThreadPool pool )
    { reinterpret_cast<WorkerMethod*>(arg)->run_p(pool); }
};

}

inline
void Run( ThreadPool pool ,
          void (*func)( void * , ThreadPool ) ,
          void * arg )
{
  TPI_Run( pool, reinterpret_cast< TPI_parallel_subprogram >(func), arg );
}

template<class Worker>
inline
void Run( ThreadPool pool, Worker & worker, void (Worker::*method)(ThreadPool))
{
  typedef WorkerMethod<Worker> WM ;

  WM tmp( worker , method );

  TPI_Run( pool, reinterpret_cast<TPI_parallel_subprogram>(& WM::run), &tmp);
}

//----------------------------------------------------------------------

}

#endif

