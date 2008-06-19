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

#include <stdlib.h>
#include <memory.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <mesh/Kernel.hpp>
#include <mesh/Entity.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Schema.hpp>
#include <mesh/FieldData.hpp>

namespace phdmesh {

namespace {

void memory_copy( unsigned char * dst , unsigned char * src , unsigned n )
{ memcpy( dst , src , n ); }

void memory_zero( unsigned char * dst , unsigned n )
{ memset( dst , 0 , n ); }

}

//----------------------------------------------------------------------
// KernelKey key = ( part-count , { part-ordinals } , counter )
//  key[ key[0] ] == counter

namespace {

unsigned kernel_counter( const unsigned * const key )
{ return key[ *key ]; }

// The part count and parts are equal
bool kernel_equal( const unsigned * lhs , const unsigned * rhs )
{
  bool result = true ;
  {
    const unsigned * const end_lhs = lhs + *lhs ;
    while ( result && end_lhs != lhs ) {
      result = *lhs == *rhs ;
      ++lhs ; ++rhs ;
    }
  }
  return result ;
}

}

// The part count and part ordinals are less
bool KernelLess::operator()( const unsigned * lhs ,
                             const unsigned * rhs ) const
{
  const unsigned * const last_lhs = lhs + ( *lhs < *rhs ? *lhs : *rhs );
  while ( last_lhs != lhs && *lhs == *rhs ) { ++lhs ; ++rhs ; }
  return *lhs < *rhs ;
}

bool Kernel::has_superset( const Part & p ) const
{
  const unsigned ordinal = p.schema_ordinal();
  const unsigned *       key_beg = key();
  const unsigned * const key_end = key_beg + *key_beg ; ++key_beg ;

  key_beg = std::lower_bound( key_beg , key_end , ordinal );

  return key_beg < key_end && ordinal == *key_beg ;
}

bool Kernel::has_superset( const PartSet & ps ) const
{
  const unsigned *       key_beg = key();
  const unsigned * const key_end = key_beg + *key_beg ; ++key_beg ;

  bool result = ! ps.empty();

  for ( PartSet::const_iterator
        i = ps.begin() ; result && i != ps.end() ; ++i ) {

    const unsigned ordinal = (*i)->schema_ordinal();

    key_beg = std::lower_bound( key_beg , key_end , ordinal );

    result = key_beg < key_end && ordinal == *key_beg ;
  }
  return result ;
}

void Kernel::supersets( PartSet & ps ) const
{
  const Schema & schema = m_mesh.schema();
  const PartSet & parts = schema.get_parts();

  const unsigned * key_val = key();
  const unsigned n = *key_val - 1 ; ++key_val ;

  ps.resize( n );
  for ( unsigned i = 0 ; i < n ; ++i , ++key_val ) {
    ps[i] = parts[ *key_val ];
  }
}

void Kernel::supersets( std::vector<unsigned> & ps ) const
{
  const unsigned * key_val = key();
  const unsigned n = *key_val - 1 ; ++key_val ;

  ps.resize( n );
  for ( unsigned i = 0 ; i < n ; ++i , ++key_val ) {
    ps[i] = *key_val ;
  }
}

//----------------------------------------------------------------------

bool field_data_valid( const FieldBase & f ,
                       const Kernel & k ,
                       unsigned ord ,
                       const char * required_by )
{
  const Schema * const k_schema = & k.mesh().schema();
  const Schema * const f_schema = & f.schema();
  const bool ok_schema  = k_schema == f_schema ;
  const bool ok_ord     = ord < k.size() ;
  const bool exists     = ok_schema && ok_ord && NULL != field_data( f , k );

  if ( required_by && ! exists ) {
    std::ostringstream msg ;
    msg << "phdmesh::field_data_valid( " ;
    msg << f ;
    msg << " , " ;
    msg << k ;
    msg << " , " ;
    msg << ord ;
    msg << " , " ;
    msg << required_by ;
    msg << " ) FAILED with " ;
    if ( ! ok_schema ) {
      msg << " different Schema" ;
    }
    else if ( ! ok_ord ) {
      msg << " Ordinal " ;
      msg << ord ;
      msg << " >= " ;
      msg << " size " ;
      msg << k.size();
    }
    else {
      msg << " no data" ;
    }
    throw std::runtime_error( msg.str() );
  }

  return exists ;
}

//----------------------------------------------------------------------

Kernel::Kernel( Mesh & arg_mesh ,
                EntityType arg_type ,
                const unsigned * arg_key )
: SetvMember<const unsigned * const>( arg_key ),
  m_mesh( arg_mesh ) ,
  m_entity_type( arg_type ),
  m_size( 0 ),
  m_capacity( 0 ),
  m_entities( NULL )
{}


//----------------------------------------------------------------------

void Kernel::zero_fields( Kernel & k_dst , unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().schema().get_fields();

  unsigned char * const p = reinterpret_cast<unsigned char*>(k_dst.m_entities);
  const DataMap *       i = k_dst.m_field_map ;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i ) {
    if ( i->m_size ) {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void Kernel::copy_fields( Kernel & k_dst , unsigned i_dst ,
                          Kernel & k_src , unsigned i_src )
{
  static const char method[] = "phdmesh::Kernel::copy_fields" ;

  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().schema().get_fields();

  unsigned char * const s = reinterpret_cast<unsigned char*>(k_src.m_entities);
  unsigned char * const d = reinterpret_cast<unsigned char*>(k_dst.m_entities);
  DataMap *       j = k_src.m_field_map ;
  DataMap *       i = k_dst.m_field_map ;
  DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i , ++j ) {

    if ( i->m_size ) {
      if ( j->m_size ) {
        if ( i->m_size == j->m_size ) {
          memory_copy( d + i->m_base + i->m_size * i_dst ,
                       s + j->m_base + j->m_size * i_src , i->m_size );
        }
        else {
          std::ostringstream msg ;
          msg << method ;
          msg << " FAILED WITH INCOMPATIBLE FIELD SIZES" ;
          throw std::runtime_error( msg.str() );
        }
      }
      else {
        memory_zero( d + i->m_base + i->m_size * i_dst , i->m_size );
      }
    }
  }
}

//----------------------------------------------------------------------

namespace {

inline size_t align( size_t nb )
{
  enum { BYTE_ALIGN = 16 };
  const size_t gap = nb % BYTE_ALIGN ;
  if ( gap ) { nb += BYTE_ALIGN - gap ; }
  return nb ;
}

const FieldBase::Dim & dimension( const FieldBase & field ,
                                  EntityType etype ,
                                  const PartSet & parts ,
                                  const char * const method )
{
  static const FieldBase::Dim empty ;

  const FieldBase::Dim * dim = & empty ;

  for ( PartSet::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const FieldBase::Dim & tmp = field.dimension( etype , **i );

    if ( tmp.part ) {
      if ( NULL == dim->part ) { dim = & tmp ; } 

      if ( Compare< MaximumFieldDimension >::
             not_equal( tmp.stride , dim->stride ) ) {

        std::ostringstream msg ;
        msg << method ;
        msg << " FAILED WITH INCOMPATIBLE DIMENSIONS FOR " ;
        msg << field ;
        msg << " Part[" << tmp.part->name() ;
        msg << "] and Part[" << dim->part->name() ;
        msg << "]" ;
     
        throw std::runtime_error( msg.str() );
      }
    }
  }

  return *dim ;
}

}

//----------------------------------------------------------------------

void Kernel::update_state()
{
  if ( 0 == kernel_counter( key() ) ) {

    const Schema & S = m_mesh.schema();
    const std::vector<FieldBase*> & field_set = S.get_fields();

    for ( unsigned i = 0 ; i < field_set.size() ; ) {

      DataMap * const tmp = m_field_map + i ;
      const FieldBase & field = * field_set[i] ;
      const unsigned num_state = field.number_of_states();
      i += num_state ;

      if ( 1 < num_state && tmp->m_size ) {
        unsigned offset[ MaximumFieldStates ] ;

        for ( unsigned j = 0 ; j < num_state ; ++j ) {
          offset[j] = tmp[j].m_base ;
        }

        for ( unsigned j = 0 ; j < num_state ; ++j ) {
          const unsigned j_new = ( j + num_state - 1 ) % num_state ;
          tmp[j_new].m_base = offset[j] ;
        }
      }
    }
  }
}

//----------------------------------------------------------------------

Kernel::~Kernel()
{}

void Mesh::destroy_kernel( KernelSet::iterator ik )
{
  Kernel * const k = & * ik ;

  KernelSet & kernel_set = m_kernels[ k->entity_type() ];

  kernel_set.remove( *ik ); // 'ik' is now invalidated

  if ( 0 == kernel_counter( k->key() ) ) {
    free( (void*) k->m_field_map );
  }
  k->m_field_map = NULL ;

  k->~Kernel();

  free( (void*) k );
}

//----------------------------------------------------------------------

KernelSet::iterator
Mesh::declare_kernel( const EntityType arg_entity_type ,
                      const PartSet & entity_parts )
{
  static const char method[] = "phdmesh::Kernel" ;
  const unsigned max = ~((unsigned) 0);

  KernelSet & kernel_set = m_kernels[ arg_entity_type ];

  const std::vector< FieldBase * > & field_set = m_schema.get_fields();

  // Entity parts are complete and contain all supersets
  // Initial key is upper bound key for kernel with these parts

  const unsigned key_size = entity_parts.size() + 2 ;

  std::vector<unsigned> key_tmp( key_size );
  key_tmp[0] = key_size - 1 ;
  key_tmp[ key_tmp[0] ] = max ;
  for ( unsigned i = 0 ; i < entity_parts.size() ; ++i ) {
    key_tmp[i+1] = entity_parts[i]->schema_ordinal();
  }

  // Does this kernel already exist?

  const unsigned * const key = & key_tmp[0] ;

  KernelSet::iterator ik = kernel_set.upper_bound( key );

  // Key value for first kernel with these parts

  key_tmp[ key[0] ] = 0 ;

  Kernel * const match_kernel =
    ( ik != kernel_set.begin() ) && kernel_equal( (--ik)->key() , key )
    ? ( & * ik ) : NULL ;

  Kernel * kernel = NULL ;

  if ( match_kernel != NULL ) {
    const size_t cap = ik->capacity();
    if ( ik->size() < cap ) {
      kernel = match_kernel ;
    }
    else if ( match_kernel->key()[ key[0] ] < max ) {
      key_tmp[ key[0] ] = 1 + match_kernel->key()[ key[0] ] ;
    }
    else {
      // ERROR insane number of kernels!
      std::string msg ;
      msg.append( method );
      msg.append( " FAILED due to impossibly large number of kernels" );
      throw std::logic_error( msg );
    }
  }

  // If it does not exist must allocate and insert

  if ( kernel == NULL ) {

    const unsigned num_fields = field_set.size();

    Kernel::DataMap * field_map = NULL ;

    if ( match_kernel != NULL ) {
      field_map = match_kernel->m_field_map ;
    }
    else {
      field_map = reinterpret_cast<Kernel::DataMap*>(
                    malloc( sizeof(Kernel::DataMap) * ( num_fields + 1 ) ) );

      size_t value_offset = 0 ;

      value_offset += align( sizeof(Entity*) * m_kernel_capacity );

      for ( unsigned i = 0 ; i < num_fields ; ++i ) {
        const FieldBase  & field = * field_set[i] ;

        unsigned value_size = 0 ;

        const FieldBase::Dim & dim =
          dimension( field, arg_entity_type , entity_parts, method);

        if ( dim.part ) { // Exists

          const unsigned scalar_size =
            NumericEnum<void>::size( field.numeric_type_ordinal() );

          const unsigned field_num_dim = field.number_of_dimensions();

          value_size = scalar_size *
            ( field_num_dim ? dim.stride[ field_num_dim - 1 ] : 1 );
        }

        field_map[i].m_base = value_offset ;
        field_map[i].m_size = value_size ;
        field_map[i].m_stride = dim.stride ;

        value_offset += align( value_size * m_kernel_capacity );
      }
      field_map[ num_fields ].m_base  = value_offset ;
      field_map[ num_fields ].m_size = 0 ;
      field_map[ num_fields ].m_stride = NULL ;
    }

    // Allocation size:
    //   sizeof(Kernel) +
    //   key_size * sizeof(unsigned) +
    //   sizeof(Entity*) * capacity() +
    //   sum[number_of_fields]( fieldsize * capacity )

    const size_t size = align( sizeof(Kernel) ) +
                        align( sizeof(unsigned) * key_size ) +
                        field_map[ num_fields ].m_base ;

    // All fields checked and sized, Ready to allocate

    unsigned char * ptr = (unsigned char *) malloc( size );

    kernel = (Kernel *) ptr ; ptr += align( sizeof(Kernel) );

    {
      unsigned * new_key = (unsigned *) ptr ;

      ptr += align( sizeof(unsigned) * key_size );

      for ( unsigned i = 0 ; i < key_size ; ++i ) { new_key[i] = key[i] ; }

      new(kernel) Kernel( *this , arg_entity_type , new_key );
    }

    kernel->m_size      = 0 ;
    kernel->m_capacity  = m_kernel_capacity ;
    kernel->m_field_map = field_map ;
    kernel->m_entities  = (Entity **) ptr ;

    std::pair<KernelSet::iterator,bool> result = kernel_set.insert( kernel );

    if ( ! result.second ) {
      std::string msg ;
      msg.append( method );
      msg.append( " FAILED INSERTION" );
      throw std::logic_error( msg );
    }

    ik = result.first ;
  }
  return ik ;
}

//----------------------------------------------------------------------

void Mesh::remove_entity( KernelSet::iterator ik , unsigned i )
{
  const unsigned entity_type = ik->m_entity_type ;

  KernelSet & kset = m_kernels[ entity_type ] ;

  // Find the last compatible kernel and move the
  // last entity up to fill in the gap.
  // Compatible entities have the same part list (key.first)
  // and incrementing ordinals (key.second).
  // Thus the next kernel with a zero ordinal will have
  // a different part list.

  const KernelSet::iterator ek = kset.end();

  KernelSet::iterator jk = ik ;
  while ( ++jk != ek && kernel_counter( jk->key() ) );
  --jk ;

  // Only move if not the last entity being removed

  if ( jk != ik || ik->m_size != i + 1 ) {

    // Not the same kernel or not the last entity

    // Copy last entity in jk to ik slot i

    Kernel::copy_fields( *ik , i , *jk , jk->m_size - 1 );

    ik->m_entities[i] = jk->m_entities[ jk->m_size - 1 ];
    ik->m_entities[i]->m_kernel     = ik ;
    ik->m_entities[i]->m_kernel_ord = i ;
  }

  --( jk->m_size );

  if ( jk->m_size != 0 ) {
    jk->m_entities[ jk->m_size ] = NULL ;
  }
  else {
    destroy_kernel( jk );
  }
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Kernel & k )
{
  const char * const entity_name = entity_type_name( k.entity_type() );

  PartSet parts ; k.supersets( parts );

  s << "Kernel( " << entity_name << " : " ;
  for ( PartSet::iterator i = parts.begin() ; i != parts.end() ; ++i ) {
    s << (*i)->name() << " " ;
  }
  s << ")" ;

  return s ;
}


std::ostream &
Kernel::print( std::ostream & os , const std::string & lead ) const
{
  const Schema & schema = m_mesh.schema();
  const PartSet & parts = schema.get_parts();
  const char * const entity_name = entity_type_name( m_entity_type );

  const unsigned * k = key();
  const unsigned n = *k ; ++k ;

  os << lead
     << "Kernel(" << entity_name << " : " ;
  for ( unsigned i = 1 ; i < n ; ++i , ++k ) {
    const std::string & name = parts[ *k ]->name(); os << " " << name ;
  }

  os << " " << *k << " ){ "
     << m_size << " of " << m_capacity << " }"
     << std::endl << lead
     << "  members {" ;

  for ( unsigned j = 0 ; j < m_size ; ++j ) {
    const long id = m_entities[j]->identifier(); os << " " << id ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

