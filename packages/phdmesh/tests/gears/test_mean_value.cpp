/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
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

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>

#include <util/TPI.h>
#include <util/SimpleArrayOps.hpp>
#include <util/ParallelComm.hpp>

#include <mesh/FieldTraits.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/FieldParallel.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/Comm.hpp>

#include <element/Hexahedron_Topologies.hpp>
#include <element/Dimensions.hpp>

using namespace phdmesh ;

namespace {

//----------------------------------------------------------------------

template< unsigned Length >
void mean_value( const size_t count ,
                 const size_t number ,
                 const double * const * const input ,
                 double * const output )
{
  const size_t total_out = count * Length ;

  for ( size_t k = 0 ; k < total_out ; ++k ) { output[k] = 0 ; }

  for ( size_t i = 0 ; i < count ; ++i ) {
    const double * const * const in  = input  + i * number ;
          double * const         out = output + i * Length ;

    for ( size_t j = 0 ; j < number ; ++j ) {
      Sum<Length>( out , in[j] );
    }
  }

  for ( size_t k = 0 ; k < total_out ; ++k ) { output[k] /= number ; }
}

struct ElementMeanValueOp {

  BulkData                     & m_mesh ;
  const Field<double,Cartesian>    & m_field ;
  const Field<double*,ElementNode> & m_field_ptr ;

  static void run( void * , TPI_ThreadPool );

  ElementMeanValueOp( const Field<double,Cartesian>     & arg_field ,
                      const Field<double*,ElementNode>  & arg_field_ptr ,
                      BulkData & arg_mesh );
};

ElementMeanValueOp::ElementMeanValueOp(
  const Field<double,Cartesian>     & arg_field ,
  const Field<double*,ElementNode>  & arg_field_ptr ,
  BulkData & arg_mesh )
: m_mesh( arg_mesh ),
  m_field( arg_field ),
  m_field_ptr( arg_field_ptr )
{
   TPI_Run( reinterpret_cast<TPI_parallel_subprogram>( & ElementMeanValueOp::run ) , this , 0 );
}

void ElementMeanValueOp::run( void * arg , TPI_ThreadPool pool )
{
  int p_rank , p_size ;

  TPI_Rank( pool , & p_rank , & p_size );

  const ElementMeanValueOp & op = * reinterpret_cast<ElementMeanValueOp*>(arg);

  Part & owns = op.m_mesh.mesh_meta_data().locally_owned_part();

  const KernelSet & ks = op.m_mesh.kernels( Element );

  KernelSet::const_iterator ik = ks.begin();

  // Advance to this thread's starting point:
  for ( int i = 0 ; ik != ks.end() && i < p_rank ; ++i , ++ik );

  while ( ik != ks.end() ) {

    if ( ik->has_superset( owns ) ) {
      const Kernel & k = *ik ;

      KernelArray< Field<double*,ElementNode> > node_array(op.m_field_ptr, k);
      KernelArray< Field<double,Cartesian> > elem_array( op.m_field , k );

      switch( elem_array.dimension<0>() ) {
      case 3 :
        mean_value<3>( node_array.dimension<1>() ,
                       node_array.dimension<0>() ,
                       node_array.contiguous_data() ,
                       elem_array.contiguous_data() );
        break ;
      }
    }

    // Advance to this thread's next kernel:
    for ( int i = 0 ; ik != ks.end() && i < p_size ; ++i , ++ik );
  }
}

}

//----------------------------------------------------------------------

namespace {

struct NodeMeanValueOp {

  BulkData                  & m_mesh ;
  const Field<double,Cartesian> & m_node_field ;
  const Field<double,Cartesian> & m_elem_field ;

  static void run( void * , TPI_ThreadPool );

  NodeMeanValueOp( const Field<double,Cartesian> & arg_node_field ,
                   const Field<double,Cartesian> & arg_elem_field ,
                   BulkData & arg_mesh );
};

NodeMeanValueOp::NodeMeanValueOp(
  const Field<double,Cartesian> & arg_node_field ,
  const Field<double,Cartesian> & arg_elem_field ,
  BulkData & arg_mesh )
: m_mesh( arg_mesh ),
  m_node_field( arg_node_field ),
  m_elem_field( arg_elem_field )
{
   TPI_Run( reinterpret_cast<TPI_parallel_subprogram>( & NodeMeanValueOp::run ) , this , 0 );
}

void NodeMeanValueOp::run( void * arg , TPI_ThreadPool pool )
{
  int p_rank , p_size ;

  TPI_Rank( pool , & p_rank , & p_size );

  const NodeMeanValueOp & op = * reinterpret_cast<NodeMeanValueOp*>(arg);

  Part & uses = op.m_mesh.mesh_meta_data().locally_used_part();

  const KernelSet & ks = op.m_mesh.kernels( Node );

  KernelSet::const_iterator ik = ks.begin();

  // Advance to this thread's starting point:
  for ( int i = 0 ; ik != ks.end() && i < p_rank ; ++i , ++ik );

  while ( ik != ks.end() ) {

    if ( ik->has_superset( uses ) ) {
      const Kernel & k = *ik ;

      KernelArray< Field<double,Cartesian> > node_array( op.m_node_field , k );

      const unsigned length = node_array.dimension<0>();
      const unsigned count  = node_array.dimension<1>();

      for ( unsigned i = 0 ; i < count ; ++i ) {
        const PairIterRelation rel = k[i]->relations( Element );
        const unsigned num_elem = rel.size();

        for ( unsigned j = 0 ; j < num_elem ; ++j ) {
          Entity & elem = * rel[j].entity();

          EntityArray< Field<double,Cartesian> >
            elem_array( op.m_elem_field , elem );

          for ( unsigned m = 0 ; m < length ; ++m ) {
            node_array(m,i) += elem_array(m);
          }
        }

        for ( unsigned m = 0 ; m < length ; ++m ) {
          node_array(m,i) /= num_elem ;
        }
      }
    }

    // Advance to this thread's next kernel:
    for ( int i = 0 ; ik != ks.end() && i < p_size ; ++i , ++ik );
  }
}

}

//----------------------------------------------------------------------
// Average a field to the element and then back to the node:

void test_diffuse_field(
  BulkData & mesh ,
  const Field<double,Cartesian> & arg_field ,
  const ElementNodePointerField & arg_field_ptr )
{
  ElementMeanValueOp( arg_field , arg_field_ptr , mesh );

  {
    const FieldBase * const f = & arg_field ;
    
    communicate_field_data( mesh , mesh.ghost_source() ,
                                   mesh.ghost_destination() ,
                                   std::vector<const FieldBase *>(1,f));
  }

  NodeMeanValueOp( arg_field , arg_field , mesh );
}
 
