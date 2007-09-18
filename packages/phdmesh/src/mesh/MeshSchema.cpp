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

void Schema::assert_not_predefined( const char * method , Part & p ) const
{
  if ( & p == & m_universal_part ||
       & p == m_owns_part     ||
       & p == m_shares_part    ||
       & p == m_aura_part ) {
    std::string msg ;
    msg.append( method )
       .append( "( " )
       .append( p.name() )
       .append( " ) FAILED Is predefined part" );
    throw std::logic_error( msg );
  }
}

//----------------------------------------------------------------------

Schema::Schema( unsigned arg_dimension , ParallelMachine pm )
  : m_commit( false ),
    m_dimension( arg_dimension ),
    m_parallel_machine( pm ),
    m_parallel_size( parallel_machine_size( pm ) ),
    m_parallel_rank( parallel_machine_rank( pm ) ),
    m_universal_part( *this , std::string( "{UNIVERSAL}" ) , PartSet() ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_aura_part( NULL )
{
  { // Declare remaining predefined parts
    const std::string owns_part_name(   "{PARALLEL_OWNS}" );
    const std::string shares_part_name( "{PARALLEL_SHARES}" );
    const std::string aura_part_name(   "{PARALLEL_AURA}" );

    m_owns_part  = & declare_part( owns_part_name );
    m_shares_part = & declare_part( shares_part_name );
    m_aura_part   = & declare_part( aura_part_name );
  }
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

namespace {

void clean_intersection( const char * const method ,
                         PartSet     & pset ,
                         std::string & name )
{
  static const char separator[] = "^" ;

  order( pset );

  PartSet::iterator i ;

  for ( i = pset.begin() ; i != pset.end() ; ) {
    // If a subset of 'i' is contained then 'i' is redundant
    if ( intersect( (*i)->subsets() , pset ) ) {
      i = pset.erase( i );
    }
    else {
      ++i ;
    }
  }

  if ( pset.size() < 2 ) {
    std::string msg ;
    msg.append(method);
    msg.append(" : FAILED, Cannot intersect fewer than two unique parts." );
    msg.append(" Input {" );
    PartSet::iterator j ;
    for ( j = pset.begin() ; j != pset.end() ; ++j ) {
      msg.append(" ");
      msg.append( (*j)->name() );
    }
    msg.append(" } Clean {" );
    for ( i = pset.begin() ; i != pset.end() ; ) {
      msg.append(" ");
      msg.append( (*i)->name() );
    }
    throw std::invalid_argument(msg);
  }

  name.assign("{");
  for ( i = pset.begin() ; i != pset.end() ; ++i ) {
    if ( i != pset.begin() ) { name.append( separator ); }
    name.append( (*i)->name() );
  }
  name.append("}");
}

}

Part & Schema::declare_part( const PartSet & pset )
{
  static const char method[] = "phdmesh::Mesh::declare_part" ;

  assert_not_committed( method );

  PartSet pset_clean( pset );

  std::string p_name ;

  clean_intersection( method , pset_clean , p_name );

  Part * p = get_part( p_name , false );

  if ( p == NULL ) { p = new Part( *this , p_name , pset_clean ); }

  if ( pset_clean != p->intersection_of() ) {
    std::string msg ;
    msg.append(method);
    msg.append(" : FAILED, Redundant incompatible intersection.");
    throw std::invalid_argument(msg);
  }

  return *p ;
}

Part & Schema::declare_part( const std::string & p_name )
{
  static const char method[] = "phdmesh::Mesh::declare_part" ;

  assert_not_committed( method );

  Part * p = get_part( p_name , false );

  if ( p == NULL ) { p = new Part( *this , p_name , PartSet() ); }

  return *p ;
}

Part * Schema::get_part( const std::string & p_name ,
                         const char * required_by ) const
{
  const PartSet & all_parts = m_universal_part.m_subsets ;

  Part * const p = find( all_parts , p_name );

  if ( required_by && NULL == p ) { // ERROR
    static const char method[] = "phdmesh::Mesh::get_part" ;
    std::string msg ;
    msg.append( method )
       .append( "( " )
       .append( p_name )
       .append( " , " )
       .append( required_by )
       .append( " ) FAILED to find part" );
    throw std::runtime_error( msg );
  }

  return p ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Commit fields and parts

namespace {

void pack( CommBuffer & b , const PartSet & pset )
{
  PartSet::const_iterator i , j ;
  for ( i = pset.begin() ; i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const CSet    & cset      = p.cset_query();
    const PartSet & subsets   = p.subsets();
    const PartSet & intersect = p.intersection_of();

    const unsigned     name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    std::string cset_text ;
    cset.print( cset_text , ":" );
    const unsigned     cset_len = cset_text.size() + 1 ;
    const char * const cset_ptr = cset_text.c_str();

    {
      const unsigned ord = p.schema_ordinal();
      b.pack<unsigned>( ord );
    }

    b.pack<unsigned>( name_len );
    b.pack<char>( name_ptr , name_len );

    b.pack<unsigned>( cset_len );
    b.pack<char>( cset_ptr , cset_len );

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
    const CSet    & cset      = p.cset_query();
    const PartSet & subsets   = p.subsets();
    const PartSet & intersect = p.intersection_of();
    const unsigned     name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    std::string cset_text ;
    cset.print( cset_text , ":" );
    const unsigned     cset_len = cset_text.size() + 1 ;
    const char * const cset_ptr = cset_text.c_str();

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
      ok = b_tmp == cset_len ;
    }
    if ( ok ) {
      b.unpack<char>( b_text , cset_len );
      ok = 0 == strcmp( cset_ptr , b_text );
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

  // Assign ordinals to fields, grouped by entity type
  // Insure that the fields' dimension maps do not have
  // incompatible subsets or supersets.

  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    std::vector< Field<void,0> *> & fm = m_fields[i] ;
    std::vector< Field<void,0> *>::iterator j ;
    for ( j = fm.begin() ; j != fm.end() ; ++j ) {
      (*j)->clean_dimension();
    }
  }

  // Parallel verification of part consistency

  {
    const bool is_root = 0 == m_parallel_rank ;

    CommBroadcast comm( m_parallel_machine , 0 );

    if ( is_root ) pack( comm.send_buffer() , m_universal_part.subsets() );

    comm.allocate_buffer();

    if ( is_root ) pack( comm.send_buffer() , m_universal_part.subsets() );

    comm.communicate();

    if ( ! unpack_verify( comm.recv_buffer() , m_universal_part.subsets() ) ) {
      std::ostringstream msg ;
      msg << "P" << m_parallel_rank ;
      msg << ": " << method ;
      msg << " : FAILED Parallel Part consistency" ;
      throw std::logic_error( msg.str() );
    }
  }

  m_commit = true ; // Cannot add or change parts or fields now
}

//----------------------------------------------------------------------

} // namespace phdmesh

