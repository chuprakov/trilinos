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

#include <strings.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <mesh/Part.hpp>
#include <mesh/Schema.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

namespace {

struct PartLess {
  bool operator()( const Part * lhs , const Part & rhs ) const
  {
    const unsigned l = lhs->schema_ordinal();
    const unsigned r = rhs.schema_ordinal();
    return l < r ;
  }

  bool operator()( const Part * lhs , const Part * rhs ) const
  {
    const unsigned l = lhs->schema_ordinal();
    const unsigned r = rhs->schema_ordinal();
    return l < r ;
  }

  bool operator()( const Part * lhs , const std::string & rhs ) const
  {
    const char * const lhs_c_str = lhs->name().c_str();
    const char * const rhs_c_str = rhs.c_str();
    return strcasecmp( lhs_c_str , rhs_c_str ) < 0 ;
  }
};

PartSet::iterator
lower_bound( PartSet & parts , const std::string & name )
{
  const PartSet::iterator e = parts.end();
  const PartSet::iterator i = parts.begin();

  return std::lower_bound( i , e , name, PartLess() );
}

PartSet::const_iterator
lower_bound( const PartSet & parts , const std::string & name )
{
  const PartSet::const_iterator e = parts.end();
  const PartSet::const_iterator i = parts.begin();

  return std::lower_bound( i , e , name, PartLess() );
}

}

//----------------------------------------------------------------------

Part * find( const PartSet & parts , const std::string & name )
{
  const PartSet::const_iterator i = lower_bound( parts , name );

  return ( i != parts.end() ) && ( name == (*i)->name() ) ? *i : (Part*) NULL ;
}

//----------------------------------------------------------------------

bool verify( const Part & p , std::string & msg )
{
  bool ok = true ;

  // Superset/subset consistency

  const Part    & universal = p.schema().universal_part();
  const PartSet & supersets = p.supersets();
  const PartSet & subsets   = p.subsets();
  const PartSet & intersection = p.intersection_of();

  std::vector<Part*>::const_iterator i , j ;

  if ( & p == & universal ) {
    if ( ! supersets.empty() ) {
      ok = false ;
      msg.append(" ");
      msg.append( p.name() );
      msg.append( " cannot have supersets ;" );
    }
  }
  else {

    // Unversial superset with symmetry

    if ( ! contain( supersets , universal ) ) {
      ok = false ;
      msg.append(" ");
      msg.append( p.name() );
      msg.append( " not-in " );
      msg.append( universal.name() );
      msg.append( " ;" );
    }

    if ( ! contain( universal.subsets() , p ) ) {
      ok = false ;
      msg.append(" ");
      msg.append( universal.name() );
      msg.append( " not-contain " );
      msg.append( p.name() );
      msg.append( " ;" );
    }

    // Superset symmetry and noncircular

    for ( i = supersets.begin() ; i != supersets.end() ; ++i ) {
      Part & s = **i ;
      if ( ! contain( s.subsets() , p ) ) {
        ok = false ;
        msg.append( " Asymmetry " );
        msg.append( p.name() );
        msg.append( " in " );
        msg.append( s.name() );
        msg.append( " ;" );
      }

      if ( contain( subsets , s ) ) {
        ok = false ;
        msg.append( " Circular " );
        msg.append( p.name() );
        msg.append( " in " );
        msg.append( s.name() );
        msg.append( " ;" );
      }
    }

    // Subset symmetry, noncircular, and transitive

    for ( i = subsets.begin() ; i != subsets.end() ; ++i ) {
      Part & sub = **i ;
      if ( ! contain( sub.supersets() , p ) ) {
        ok = false ;
        msg.append( " Asymmetry " );
        msg.append( p.name() );
        msg.append( " contain " );
        msg.append( sub.name() );
        msg.append( " ;" );
      }
      if ( contain( supersets , sub ) ) {
        ok = false ;
        msg.append( " Circular " );
        msg.append( p.name() );
        msg.append( " contain " );
        msg.append( sub.name() );
        msg.append( " ;" );
      }
      for ( j = supersets.begin() ; j != supersets.end() ; ++j ) {
        Part & s = **j ;
        if ( ! contain( s.subsets() , sub ) ) {
          ok = false ;
          msg.append( " Not-transitive " );
          msg.append( s.name() );
          msg.append( " contain " );
          msg.append( p.name() );
          msg.append( " contain " );
          msg.append( sub.name() );
          msg.append( " ;" );
        }
      }
    }
  }

  // Intersection is subclass of supersets.
  // Intersection members cannot be subset/superset of one another.

  for ( i = intersection.begin() ; i != intersection.end() ; ) {
    Part & a = **i ; ++i ;
    if ( ! contain( supersets , a ) ) {
      ok = false ;
      msg.append( " Intersection-superset " );
      msg.append( p.name() );
      msg.append( " not-in " );
      msg.append( a.name() );
      msg.append( " ;" );
    }

    for ( j = i ; j != intersection.end() ; ) {
      Part & b = **j ; ++j ;
      if ( contain( a.supersets() , b ) ||
           contain( b.supersets() , a ) ) {
        ok = false ;
        msg.append( " Intersection " );
        msg.append( p.name() );
        msg.append( " ( " );
        msg.append( a.name() );
        msg.append( " , " );
        msg.append( b.name() );
        msg.append( " );" );
      }
    }
  }

  return ok ;
}

//----------------------------------------------------------------------

std::ostream &
print( std::ostream & os , const char * const lead , const Part & p )
{
  const PartSet & supersets = p.supersets();
  const PartSet & subsets   = p.subsets();
  const PartSet & intersection = p.intersection_of();

  std::vector<Part*>::const_iterator i ;

  if ( lead != NULL ) { os << lead ; }
  os << "Part[ " ;
  os << p.name() ;
  os << " , " ;
  os << p.schema_ordinal() ;
  os << " ] {" ;
  os << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Supersets {" ;
  for ( i = supersets.begin() ; i != supersets.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " }" << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Intersection_Of {" ;
  for ( i = intersection.begin() ; i != intersection.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " } }" << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Subsets {" ;
  for ( i = subsets.begin() ; i != subsets.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

void order( PartSet & v )
{
  PartSet::iterator ev = v.end();
  PartSet::iterator iv = v.begin();
  std::sort( iv , ev , PartLess() );
  iv = std::unique( iv , ev );
  v.erase( iv , ev );
}

bool insert( PartSet & v , Part & part )
{
  const PartSet::iterator e = v.end();
        PartSet::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  const bool new_member = i == e || *i != & part ;

  if ( new_member ) { Part * const tmp = & part ; v.insert( i , tmp ); }
  return new_member ;
}

bool contain( const PartSet & v , const Part & part )
{
  const PartSet::const_iterator e = v.end();
        PartSet::const_iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  return i != e && *i == & part ;
}

bool contain( const PartSet & super , const PartSet & sub )
{
  bool result = ( ! sub.empty() ) && ( sub.size() <= super.size() );

  if ( result ) {
    PartLess comp ;

    const PartSet::const_iterator ev = super.end();
          PartSet::const_iterator iv = super.begin();

    const PartSet::const_iterator ep = sub.end();
          PartSet::const_iterator ip = sub.begin();

    while ( result && ip != ep ) {
      Part * const q = *ip ; ++ip ;
      iv = std::lower_bound( iv , ev , q , comp );
      result = iv != ev && *iv == q ; 
    }
  }

  return result ;
}

unsigned intersect( const PartSet & v , const PartSet & p )
{
  // Both lists must be sorted, assume v.size() > p.size()

  const PartSet::const_iterator ev = v.end();
        PartSet::const_iterator iv = v.begin();

  const PartSet::const_iterator ep = p.end();
        PartSet::const_iterator ip = p.begin();

  unsigned count = 0 ;

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { ++count ; }
  }

  return count ;
}

unsigned intersect( const PartSet & v , const PartSet & p , PartSet & r )
{
  // Both lists must be sorted, assume v.size() > p.size()

  const PartSet::const_iterator ev = v.end();
        PartSet::const_iterator iv = v.begin();

  const PartSet::const_iterator ep = p.end();
        PartSet::const_iterator ip = p.begin();

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { r.push_back( q ); }
  }

  return r.size() ;
}

bool intersect( const Part & a , const Part & b )
{
  const PartSet & a_sub = a.subsets();
  const PartSet & b_sub = b.subsets();
  return contain( a_sub , b ) ||
         contain( b_sub , a ) ||
         intersect( b_sub , a_sub );
}

//----------------------------------------------------------------------

Part::~Part()
{}

Part::Part( Schema & m , const std::string & n , const PartSet & intersect )
  : m_name( n ),
    m_cset(),
    m_schema( m ) ,
    m_schema_ordinal( 0 ) ,
    m_subsets() , m_supersets() , m_intersect()
{
  static const char method[] = "phdmesh::Part::Part" ;

  Part * const universal = const_cast<Part*>( & m.universal_part() );

  if ( universal != this ) {
    m_supersets.push_back( universal );
  }

  PartSet & all_parts = universal->m_subsets ;

  {
    PartSet::iterator i = lower_bound( all_parts , n );

    if ( i != all_parts.end() && (*i)->m_name == n ) {
      std::string msg ;
      msg.append( method );
      msg.append( "( <schema> , " );
      msg.append( n );
      msg.append( " ) : FAILED, duplicate part in schema." );
      throw std::invalid_argument( msg );
    }

    Part * const tmp = this ;
    all_parts.insert( i , tmp );
  }

  // Update the mesh parts' ordinals
  for ( unsigned j = 0 ; j < all_parts.size() ; ++j ) {
    all_parts[j]->m_schema_ordinal = j ;
  }

  if ( ! intersect.empty() ) {
    m_intersect = intersect ;
    PartSet::const_iterator i ;
    for ( i = intersect.begin() ; i != intersect.end() ; ++i ) {
      m_schema.declare_part_subset( **i , *this );
    }
  }
}

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

  Part * p = find( m_universal_part.m_subsets , p_name );

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

  Part * p = find( m_universal_part.m_subsets , p_name );

  if ( p == NULL ) { p = new Part( *this , p_name , PartSet() ); }

  return *p ;
}

//----------------------------------------------------------------------

void Schema::declare_part_subset( Part & superset , Part & subset )
{
  static const char method[] = "phdmesh::Schema::declare_part_subset" ;

  if ( ! contain( superset.m_subsets , subset ) ) {

    assert_not_committed( method );
    assert_same_schema(   method , superset.schema() );
    assert_same_schema(   method , subset.schema() );

    if ( & m_universal_part == & subset ||
         & superset         == & subset ||
         contain( superset.m_supersets , subset ) ) {
      std::string msg ;
      msg.append( method )
         .append( "[ " )
         .append( superset.name() )
         .append( " ] ( " )
         .append( subset.name() )
         .append( " ) FAILED, IS CIRCULAR" );
      throw std::invalid_argument( msg );
    }

    // Symmetry:

    insert( subset.m_supersets , superset );
    insert( superset.m_subsets , subset );

    PartSet::iterator i ;

    // Transitive, is also a subset of superset's supersets:

    for ( i =  superset.m_supersets.begin() ;
          i != superset.m_supersets.end() ; ++i ) {
      declare_part_subset( **i , subset );
    }

    // Intersection-part membership check.

    for ( i =  superset.m_subsets.begin() ;
          i != superset.m_subsets.end() ; ++i ) {
      Part & p_sub = **i ;
      if ( & p_sub != & subset ) {

        // If subset is fully contained in a superset->subset
        // intersection-part then is also a subset of that
        // intersection-part.

        // If subset is an intersection-part and one of the superset's
        // subsets is fully contained in the subset's intersection then
        // that superset-subset is a subset of this subset.

        if ( contain( subset.m_supersets , p_sub.m_intersect ) ) {
          declare_part_subset( p_sub , subset );
        }
        else if ( contain( p_sub.m_supersets , subset.m_intersect ) ) {
          declare_part_subset( subset , p_sub );
        }
      }
    }
  }
}

//----------------------------------------------------------------------

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

} // namespace phdmesh

