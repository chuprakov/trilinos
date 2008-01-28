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

#include <string.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <mesh/Schema.hpp>
#include <mesh/Comm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

void Schema::assert_not_committed( const char * method ) const
{
  if ( m_commit ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: mesh Schema has been committed." );
    throw std::logic_error( msg );
  }
}

void Schema::assert_committed( const char * method ) const
{
  if ( ! m_commit ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: mesh Schema has not been committed." );
    throw std::logic_error( msg );
  }
}

void Schema::assert_same_schema( const char * method ,
                                 const Schema & rhs ) const
{
  if ( this != & rhs ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED Different schema." );
    throw std::logic_error( msg );
  }
}

//----------------------------------------------------------------------

Schema::Schema()
  : m_commit( false ),
    m_universal_part( *this , std::string( "{UNIVERSAL}" ) , PartSet() ),
    m_uses_part( NULL ),
    m_owns_part( NULL )
{
  // Declare remaining predefined parts
  const std::string uses_part_name( "{USES}" );
  const std::string owns_part_name( "{OWNS}" );

  m_uses_part = & declare_part( uses_part_name );
  m_owns_part = & declare_part( owns_part_name );

  declare_part_subset( * m_uses_part , * m_owns_part );
}

void Schema::commit()
{
  static const char method[] = "phdmesh::Schema::commit" ;

  assert_not_committed( method );

  { // Verify parts
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED : " );
    const std::vector<Part*> & parts = m_universal_part.subsets();
    for ( unsigned i = 0 ; i < parts.size() ; ++i ) {
      Part & p = * parts[i] ;
      if ( ! verify( p , msg ) ) {
        throw std::logic_error( msg );
      }
    }
  }

  clean_field_dimension();

  m_commit = true ; // Cannot add or change parts or fields now
}

Schema::~Schema()
{
  // Destroy the fields, used 'new' to allocate so now use 'delete'

  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    std::vector<Field<void,0> * > & fset = m_fields[i] ;

    std::vector<Field<void,0> * >::iterator j = fset.begin();

    for ( ; j != fset.end() ; ++j ) {
      delete *j ;
    }

    fset.clear();
  }

  // Destroy the parts, used 'new' to allocate so now use 'delete'
  {
    std::vector<Part*> & parts = m_universal_part.m_subsets ;

    std::vector< Part * >::iterator j = parts.begin();

    for ( ; j != parts.end() ; ++j ) {
      if ( *j != & m_universal_part ) {
        delete *j ;
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Verify parallel consistency of fields and parts

namespace {

void pack( CommBuffer & b , const PartSet & pset )
{
  PartSet::const_iterator i , j ;
  for ( i = pset.begin() ; i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const PartSet & subsets   = p.subsets();
    const PartSet & intersect = p.intersection_of();

    const unsigned     name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    {
      const unsigned ord = p.schema_ordinal();
      b.pack<unsigned>( ord );
    }

    b.pack<unsigned>( name_len );
    b.pack<char>( name_ptr , name_len );

    const unsigned subset_size = subsets.size();
    b.pack<unsigned>( subset_size );
    for ( j = subsets.begin() ; j != subsets.end() ; ++j ) {
      const Part & s = **j ;
      const unsigned ord = s.schema_ordinal();
      b.pack<unsigned>( ord );
    }
    const unsigned intersect_size = intersect.size();
    b.pack<unsigned>( intersect_size );
    for ( j = intersect.begin() ; j != intersect.end() ; ++j ) {
      const Part & s = **j ;
      const unsigned ord = s.schema_ordinal();
      b.pack<unsigned>( ord );
    }
  }
}

bool unpack_verify( CommBuffer & b , const PartSet & pset )
{
  enum { MAX_TEXT_LEN = 4096 };
  char b_text[ MAX_TEXT_LEN ];
  unsigned b_tmp ;

  bool ok = true ;
  PartSet::const_iterator i , j ;
  for ( i = pset.begin() ; ok && i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const PartSet & subsets   = p.subsets();
    const PartSet & intersect = p.intersection_of();
    const unsigned     name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == p.schema_ordinal();
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == name_len ;
    }
    if ( ok ) {
      b.unpack<char>( b_text , name_len );
      ok = 0 == strcmp( name_ptr , b_text );
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == subsets.size() ;
    }
    for ( j = subsets.begin() ; ok && j != subsets.end() ; ++j ) {
      const Part & s = **j ;
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == s.schema_ordinal();
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == intersect.size();
    }
    for ( j = intersect.begin() ; ok && j != intersect.end() ; ++j ) {
      const Part & s = **j ;
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == s.schema_ordinal();
    }
  }
  return ok ;
}

void pack( CommBuffer & ,
           const std::vector< Field<void,0> * > & )
{
}

bool unpack_verify( CommBuffer & ,
                    const std::vector< Field<void,0> * > & )
{
  bool ok = true ;
  return ok ;
}

}

//----------------------------------------------------------------------

void verify_parallel_consistency( const Schema & s , ParallelMachine pm )
{
  static const char method[] = "phdmesh::verify_parallel_consistency(Schema)" ;

  const unsigned p_rank = parallel_machine_rank( pm );

  const bool is_root = 0 == p_rank ;

  CommBroadcast comm( pm , 0 );

  if ( is_root ) {
    pack( comm.send_buffer() , s.get_parts() );
    pack( comm.send_buffer() , s.get_fields( EntityType(0) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(1) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(2) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(3) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(4) ) );
  }

  comm.allocate_buffer();

  if ( is_root ) {
    pack( comm.send_buffer() , s.get_parts() );
    pack( comm.send_buffer() , s.get_fields( EntityType(0) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(1) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(2) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(3) ) );
    pack( comm.send_buffer() , s.get_fields( EntityType(4) ) );
  }

  comm.communicate();

  int ok[6] ;

  ok[5] = unpack_verify( comm.recv_buffer() , s.get_parts() );
  ok[0] = unpack_verify( comm.recv_buffer() , s.get_fields( EntityType(0) ) );
  ok[1] = unpack_verify( comm.recv_buffer() , s.get_fields( EntityType(1) ) );
  ok[2] = unpack_verify( comm.recv_buffer() , s.get_fields( EntityType(2) ) );
  ok[3] = unpack_verify( comm.recv_buffer() , s.get_fields( EntityType(3) ) );
  ok[4] = unpack_verify( comm.recv_buffer() , s.get_fields( EntityType(4) ) );

  all_reduce( pm , ReduceMin<6>( ok ) );

  if ( ! ( ok[0] && ok[1] && ok[2] && ok[3] && ok[4] && ok[5] ) ) {
    std::ostringstream msg ;
    msg << "P" << p_rank ;
    msg << ": " << method ;
    msg << " : FAILED for:" ;
    if ( ! ok[5] ) { msg << " Parts ," ; }
    if ( ! ok[0] ) { msg << " Node Fields ," ; }
    if ( ! ok[1] ) { msg << " Edge Fields ," ; }
    if ( ! ok[2] ) { msg << " Face Fields ," ; }
    if ( ! ok[3] ) { msg << " Element Fields ," ; }
    if ( ! ok[4] ) { msg << " Other Fields ," ; }
    throw std::logic_error( msg.str() );
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

