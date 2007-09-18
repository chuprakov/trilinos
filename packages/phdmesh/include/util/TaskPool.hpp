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

#ifndef util_TaskPool_hpp
#define util_TaskPool_hpp

#include <util/taskpool.h>

namespace phdmesh {
namespace taskpool {

inline
int resize( unsigned num_tasks )
{ return phdmesh_taskpool_resize( num_tasks ); }

/** Run 'void Worker::method()'
 *  This method uses 'num_locks' locks to control access of shared data.
 *  Locks are bracketed in a block within this method.
 *    { taskpool::lock local_lock( i );
 *    } // 'local_lock' releases at destruction.
 *  Returns number of tasks actually used for success.
 *  Returns negative number if Worker::operator() throws an exception.
 */
template<class Worker>
inline
int run( Worker & , void (Worker::*)(unsigned,unsigned), unsigned num_locks );

class lock {
private:
  lock();
  lock( const lock & );
  lock & operator = ( const lock & );
  const unsigned m_value ;
  void throw_lock();
  void throw_unlock();
public:
  explicit lock( unsigned i_lock , const char * const error_info = NULL )
    : m_value( i_lock ) { phdmesh_taskpool_lock( i_lock , error_info ); }

  ~lock() { phdmesh_taskpool_unlock( m_value , NULL ); }
};

//----------------------------------------------------------------------

namespace {

template<class Worker>
class WorkerMethod {
private:
  WorkerMethod();
  WorkerMethod( const WorkerMethod & );
  WorkerMethod & operator = ( const WorkerMethod & );

  typedef void (Worker::*Method)( unsigned p_size , unsigned p_rank );
  Worker & worker ;
  Method   method ;
public:
  WorkerMethod( Worker & w , Method m ) : worker(w), method(m) {}
  void run( unsigned p_size , unsigned p_rank )
    { (worker.*method)( p_size , p_rank ); }
};

template<class Worker>
int worker_method_run( void * arg , unsigned p_size , unsigned p_rank )
{
  int result = 0 ;
  try {
    reinterpret_cast<WorkerMethod<Worker>*>(arg)->run(p_size,p_rank);
  } catch( ... ) { result = -1 ; }
  return result ;
}

}

inline
int run( int (*routine)(void *,unsigned,unsigned) ,
         void * routine_data , unsigned num_locks )
{
  return phdmesh_taskpool_run(
    reinterpret_cast< phdmesh_taskpool_routine >( routine ) ,
    routine_data , num_locks );
}

template<class Worker>
inline
int run( Worker & worker,
         void (Worker::*method)(unsigned,unsigned),
         unsigned num_locks )
{
  WorkerMethod<Worker> tmp( worker , method );

  return run( & worker_method_run<Worker>, & tmp , num_locks );
}

//----------------------------------------------------------------------

}
}

#endif

