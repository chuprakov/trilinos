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

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include <util/CSet.hpp>

namespace phdmesh {

const char * CSetMember<void>::name() const
{
  static const char n[] = "phdmesh::CSetMember<undetermined>" ;
  return n ;
}

void CSet::default_destroy( const CSetMember<void> * m )
{
  delete const_cast<CSetMember<void>*>( m );
}

void CSet::print( std::string & s , const char * separator ) const
{
  const std::vector<CSet::ValueType>::const_iterator b = m_inst.begin();
  const std::vector<CSet::ValueType>::const_iterator e = m_inst.end();
  std::vector<CSet::ValueType>::const_iterator i ;

  for ( i = b ; i != e ; ++i ) {
    if ( separator && i != b ) { s.append( separator ); }
    s.append( "( typeid(" );
    s.append( i->first->m_cset_type.name() );
    s.append( ") , " );
    s.append( i->first->name() );
    s.append( " )" );
  }
}

//----------------------------------------------------------------------

namespace {

// Comparison for sorted vector

typedef std::pair< const CSetMember<void> * ,
                   void (*)( const CSetMember<void> * ) > VT ;

struct less_cset {
  bool operator()( const VT & lhs , const std::type_info & rhs ) const ;
  bool operator()( const std::type_info & lhs , const VT & rhs ) const ;
};

// On some systems, namely AIX, std::type_info::before(...)
// has a bug where it returns true instead of false for equality.
// Thus we pay a small price on all systems to specifically
// test for and eliminate equality.

bool less_cset::operator()( const VT & lhs , const std::type_info & rhs ) const
{
  return lhs.first->m_cset_type.before( rhs ) &&
         lhs.first->m_cset_type != rhs ;
}

bool less_cset::operator()( const std::type_info & lhs , const VT & rhs ) const
{
  return lhs.before( rhs.first->m_cset_type ) &&
         lhs != rhs.first->m_cset_type ;
}

}

std::vector<CSet::ValueType>::const_iterator
CSet::m_upper_bound( std::vector<CSet::ValueType>::const_iterator i ,
                     std::vector<CSet::ValueType>::const_iterator j ,
                     const std::type_info & t )
{ return std::upper_bound( i , j , t , less_cset() ); }

std::vector<CSet::ValueType>::const_iterator
CSet::m_lower_bound( std::vector<CSet::ValueType>::const_iterator i ,
                     std::vector<CSet::ValueType>::const_iterator j ,
                     const std::type_info & t )
{ return std::lower_bound( i , j , t , less_cset() ); }

std::vector<CSet::ValueType>::iterator
CSet::m_lower_bound( std::vector<CSet::ValueType>::iterator i ,
                     std::vector<CSet::ValueType>::iterator j ,
                     const std::type_info & t )
{ return std::lower_bound( i , j , t , less_cset() ); }


//----------------------------------------------------------------------

std::pair< std::vector<VT>::const_iterator ,
           std::vector<VT>::const_iterator >
CSet::m_bounds( const std::type_info & t ) const
{
  std::vector<CSet::ValueType>::const_iterator i = m_inst.begin();
  std::vector<CSet::ValueType>::const_iterator j = m_inst.end();

  i = m_lower_bound( i , j , t );
  j = m_upper_bound( i , j , t );

  return BoundsType( i , j );
}

//----------------------------------------------------------------------

namespace {

void throw_null_value( const char * method_name ,
                       const char * type_name )
{
  std::string msg ;
  msg.append( method_name )
     .append( "<" )
     .append( type_name )
     .append( "> FAILED: Given {null} value" );
  throw std::invalid_argument( msg );
}

void throw_redundant_value( const char * method_name ,
                            const char * type_name ,
                            const char * value_name )
{
  std::string msg ;
  msg.append( method_name )
     .append( "<" )
     .append( type_name )
     .append( ">( " )
     .append( value_name )
     .append( " ) FAILED: Redundant value has a conflict" );
  throw std::invalid_argument( msg );
}

void throw_not_value( const char * method_name ,
                      const char * type_name ,
                      const char * value_name )
{
  std::string msg ;
  msg.append( method_name )
     .append( "<" )
     .append( type_name )
     .append( ">( " )
     .append( value_name )
     .append( " ) FAILED: Not currently a value" );
  throw std::invalid_argument( msg );
}

}

//----------------------------------------------------------------------

unsigned
CSet::m_insert( const CSetMember<void> * inst ,
                CSet::MemberDestroy destroy ,
                const char * const type_name )
{
  static const char method_name[] = "phdmesh::CSet::insert" ;

  if ( inst == NULL ) { throw_null_value( method_name , type_name ); }

  // Find insertion point, make sure this insert is not imcompatible

  const std::vector<CSet::ValueType>::iterator e = m_inst.end();
        std::vector<CSet::ValueType>::iterator i = m_inst.begin();
        std::vector<CSet::ValueType>::iterator j ;

  i = m_lower_bound( i , e , inst->m_cset_type );

  // Check for a duplication of this member

  bool match = false ;

  for ( j = i ; j != e &&
                j->first->m_cset_type == inst->m_cset_type &&
                ! ( match = j->first == inst ) ; ++j );

  const unsigned ordinal = j - i ;

  if ( match ) {
    if ( j->second != destroy ) {
      throw_redundant_value( method_name , type_name , inst->name() );
    }
  }
  else {
    m_inst.insert( j , CSet::ValueType( inst , destroy ) );
  }

  return ordinal ;
}

//----------------------------------------------------------------------

unsigned
CSet::m_replace( const CSetMember<void> * c_old ,
                 const CSetMember<void> * c_new ,
                 CSet::MemberDestroy destroy ,
                 const char * const type_name )
{
  static const char method_name[] = "phdmesh::CSet::replace" ;

  if ( c_old == NULL || c_new == NULL ) {
    throw_null_value( method_name , type_name );
  }

  // Search for the old and new components
  // Old component must exist.
  // New component either must not exist or be the old component

  const std::vector<CSet::ValueType>::iterator e = m_inst.end();
        std::vector<CSet::ValueType>::iterator i = m_inst.begin();
        std::vector<CSet::ValueType>::iterator j , k ;

  i = m_lower_bound( i, e, c_old->m_cset_type );

  bool old_match = false ;

  for ( j = i ; j != e &&
                j->first->m_cset_type == c_old->m_cset_type &&
                ! ( old_match = j->first == c_old ) ; ++j );

  if ( ! old_match ) {
    throw_not_value( method_name , type_name , c_old->name() );
  }

  const unsigned ordinal = j - i ;

  if ( c_old != c_new ) {

    bool new_match = false ;

    for ( k = i ; k != e &&
                  k->first->m_cset_type == c_new->m_cset_type &&
                  ! ( new_match = k->first == c_new ) ; ++k );

    if ( new_match ) {
      throw_redundant_value( method_name , type_name , c_new->name() );
    }
  }

  j->first  = c_new ;
  j->second = destroy ;

  return ordinal ;
}

//----------------------------------------------------------------------

CSet::~CSet()
{
  const std::vector<CSet::ValueType>::iterator e = m_inst.end();
        std::vector<CSet::ValueType>::iterator i = m_inst.begin();

  while ( i != e ) {
    if ( i->first && i->second ) {
      MemberDestroy d = i->second ;
      (*d)( i->first );
    }
    i->first   = NULL ;
    i->second  = NULL ;
    ++i ;
  }
}

CSet::CSet() : m_inst() {}

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------


#ifdef UNIT_TEST

namespace phdmesh {
namespace unit_test {

void DoNotDelete( const phdmesh::CSetMember<void> * m )
{
  std::cout << "DoNotDelete(" << m->name() << ")" << std::endl ;
}

void DoDelete( const phdmesh::CSetMember<void> * m )
{
  std::cout << "DoDelete(" << m->name() << ")" << std::endl ;
  delete m ;
}

class A : public phdmesh::CSetMember<A> {
public:
  virtual ~A();
};

class B : public phdmesh::CSetMember<B> {
public:
  virtual ~B();
};

class U : public A {
public:
  const char * name() const ;
};

class V : public B {
public:
  const char * name() const ;
};

class W : public B {
public:
  const char * name() const ;
};

class X : public A , public B {
public:
  const char * name() const ;
};

class Y : public A , public B {
public:
  const char * name() const ;
};

class Z {
public:
  const char * name() const ;
};

//----------------------------------------------------------------------

int cset()
{
  CSet s ;

  U * u = new U();
  V * v = new V();
  W * w = new W();
  X * x = new X();
  Y * y = new Y();

  std::cout << "s.insert<A>(u) = " << s.insert<A>(u) << std::endl ;
  std::cout << "s.insert<B>(v) = " << s.insert<B>(v) << std::endl ;
  std::cout << "s.insert<B>(w) = " << s.insert<B>(w) << std::endl ;
  std::cout << "s.insert<A>(x) = " << s.insert<A>(x) << std::endl ;
  std::cout << "s.insert<B>(x,NULL) = " << s.insert<B>(x,NULL) << std::endl ;

  std::cout
    << "s.count<A>() == " << s.count<A>() << std::endl 
    << "s.count<B>() == " << s.count<B>() << std::endl ;
  s.print( std::cout , "  " ) << std::endl ;

  std::cout << "s.replace<A>(x,x,DoNotDelete) = "
            <<  s.replace<A>(x,x,DoNotDelete) << std::endl ;
  std::cout << "s.replace<B>(x,x,DoDelete) = "
            <<  s.replace<B>(x,x,DoDelete) << std::endl ;

  std::cout
    << "s.count<A>() == " << s.count<A>() << std::endl 
    << "s.count<B>() == " << s.count<B>() << std::endl ;
  s.print( std::cout , "  " ) << std::endl ;

  std::cout << "s.replace<A>(x,y,DoNotDelete) = "
            <<  s.replace<A>(x,y,DoNotDelete) << std::endl ;
  std::cout << "s.replace<B>(x,y,DoDelete) = "
            <<  s.replace<B>(x,y,DoDelete) << std::endl ;

  delete x ; x = NULL ;

  std::cout
    << "s.count<A>() == " << s.count<A>() << std::endl 
    << "s.count<B>() == " << s.count<B>() << std::endl ;
  s.print( std::cout , "  " ) << std::endl ;

  std::cout << "s.get<A>()->name()  == " << s.get<A>()->name() << std::endl ;
  std::cout << "s.get<A>(1)->name() == " << s.get<A>(1)->name() << std::endl ;
  std::cout << "s.get<B>()->name()  == " << s.get<B>()->name() << std::endl ;
  std::cout << "s.get<B>(1)->name() == " << s.get<B>(1)->name() << std::endl ;
  std::cout << "s.get<B>(2)->name() == " << s.get<B>(2)->name() << std::endl ;

  std::vector<const B *> b_all = s.get_all<B>();
  std::cout << "b_all.size() == " << b_all.size() << std::endl ;
  for ( std::vector<const B *>::const_iterator i = b_all.begin() ;
        i != b_all.end() ; ++i ) {
    std::cout << "b_all has '" << (*i)->name()  << "'" << std::endl ;
  }

#if 0 /* Compile errors */

  Z * z = new Z();

  s.insert<A>(v);

  s.insert<Z>(z);

  s.get<Z>();

#endif

  return 0 ;
}

//----------------------------------------------------------------------

A::~A() {}
B::~B() {}

const char * U::name() const
{
  static const char n[] = "U" ;
  return n ;
}

const char * V::name() const
{
  static const char n[] = "V" ;
  return n ;
}

const char * W::name() const
{
  static const char n[] = "W" ;
  return n ;
}

const char * X::name() const
{
  static const char n[] = "X" ;
  return n ;
}

const char * Y::name() const
{
  static const char n[] = "Y" ;
  return n ;
}

const char * Z::name() const
{
  static const char n[] = "Z" ;
  return n ;
}

}
}

int main()
{
  return phdmesh::unit_test::cset();
}

#endif



