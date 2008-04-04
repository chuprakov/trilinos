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

#include <stddef.h>
#include <util/TPI.hpp>
#include <util/Basics.hpp>
#include <mesh/Assemble.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>
#include <mesh/Schema.hpp>
#include <mesh/Kernel.hpp>
#include <mesh/Entity.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

namespace {

class AssembleTask {
private:
  AssembleTask();
  AssembleTask( const AssembleTask & );
  AssembleTask & operator = ( const AssembleTask & );
public:
  const Mesh                  & mesh ;
  const unsigned                entity_type ;
  const std::vector<Assemble> & ops ;

  KernelSet::const_iterator ik_work ;

  AssembleTask( const Mesh & arg_mesh ,
                unsigned     arg_entity_type ,
                const std::vector<Assemble> & arg_ops );

  void work(TPI::ThreadPool);
};

AssembleTask::AssembleTask( const Mesh & arg_mesh ,
                            unsigned     arg_entity_type ,
                            const std::vector<Assemble> & arg_ops )
  : mesh( arg_mesh ),
    entity_type( arg_entity_type ),
    ops( arg_ops ),
    ik_work( arg_mesh.kernels( arg_entity_type ).begin() )
{
  TPI::Set_lock_size( 1 );
  TPI::Run( *this , & AssembleTask::work );
}

void AssembleTask::work(TPI::ThreadPool pool)
{
  const std::vector<Assemble>::const_iterator ia_beg = ops.begin();
  const std::vector<Assemble>::const_iterator ia_end = ops.end();
  
  const KernelSet::const_iterator ik_end = mesh.kernels( entity_type ).end();

  const Schema & schema = mesh.schema();
  Part & uses = schema.uses_part();

  for(;;) {
    // Get work:

    const Kernel * k = NULL ;
    {
      TPI::LockGuard get_work_lock(pool,0);
      for ( ; ik_work != ik_end && NULL == k ; ++ik_work ) {
        if ( ik_work->has_superset( uses ) ) {
          k = & *ik_work ;
        }
      }
    }
    if ( NULL == k ) break ;
 
    //------------------------------
    // Do work on this kernel:

    for ( std::vector<Assemble>::const_iterator
          ia = ia_beg ; ia_end != ia ; ++ia ) {
      const Assemble & assemble = *ia ;

      if ( entity_type == assemble.dst_field().entity_type() ) {

        const Field<void,0> & dst_field = assemble.dst_field();

        const unsigned data_size = k->data_size( dst_field );

        if ( data_size ) { // Exists on this kernel

          const Field<void,0> & src_field = assemble.src_field();
          const unsigned        src_type  = src_field.entity_type();

          unsigned char * dst_ptr = (unsigned char *) k->data( dst_field );

          const Kernel::iterator ie_end = k->end();
                Kernel::iterator ie     = k->begin();

          for ( ; ie_end != ie ; ++ie , dst_ptr += data_size ) {

            for ( ConnectSpan con = (*ie)->connections( src_type );
                  con ; ++con ) {
              const unsigned src_id  = con->identifier();
              void * const   src_ptr = con->entity()->data(src_field);
              assemble( dst_ptr , src_ptr , src_id );
            }
          }
        }
      }
    }
  }
}

}

void assemble( const Mesh & M , const std::vector<Assemble> ops )
{
  static const char method[] = "phdmesh::assemble" ;

  // Visit each destination/update kernel exactly once.

  const std::vector<Assemble>::const_iterator ia_end = ops.end();
  const std::vector<Assemble>::const_iterator ia_beg = ops.begin();
        std::vector<Assemble>::const_iterator ia ;

  unsigned dst_type_present[ end_entity_rank ];

  for ( unsigned t = 0 ; t < end_entity_rank ; ++t ) {
    dst_type_present[t] = 0 ;
  }

  // Parallel copy for source data

  {
    std::vector< const Field<void,0> * > update_fields ;
    update_fields.reserve( ops.size() );

    for ( ia = ia_beg ; ia_end != ia ; ) {
      const Assemble & op = *ia ; ++ia ;
      const Field<void,0> * const src_field = & op.src_field();
      const Field<void,0> * const dst_field = & op.dst_field();

      M.schema().assert_same_schema( method , dst_field->schema() );
      M.schema().assert_same_schema( method , src_field->schema() );

      dst_type_present[ dst_field->entity_type() ] = 1 ;

      update_fields.push_back( src_field );
    }

    comm_mesh_field_values( M , M.aura_domain() , M.aura_range() ,
                            update_fields , false );
  }

  // Local assembly

  for ( unsigned t = 0 ; t < end_entity_rank ; ++t ) {
    if ( dst_type_present[ t ] ) {
      AssembleTask work( M , t , ops );
    }
  }
}

}

