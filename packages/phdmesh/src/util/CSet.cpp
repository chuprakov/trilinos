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
#include <sstream>

#include <util/CSet.hpp>

namespace phdmesh {

namespace {

typedef void (* DeleteFunction )( void * );

typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

// Comparison for sorted vector

struct less_cset {
  bool operator()( const Manager        & lhs ,
                   const std::type_info & rhs ) const ;
  bool operator()( const std::type_info & lhs ,
                   const Manager        & rhs ) const ;
};

// On some systems, namely AIX, std::type_info::before(...)
// has a bug where it returns true instead of false for equality.
// Thus we pay a small price on all systems to specifically
// test for and eliminate equality.

bool less_cset::operator()( const Manager        & lhs ,
                            const std::type_info & rhs ) const
{ return lhs.first->before( rhs ) && * lhs.first != rhs ; }

bool less_cset::operator()( const std::type_info & lhs ,
                            const Manager        & rhs ) const
{ return lhs.before( *rhs.first ) && lhs != *rhs.first ; }

std::pair<size_t,size_t>
span( const std::vector< Manager > & v , const std::type_info & t )
{
  std::vector< Manager >::const_iterator i = v.begin();
  std::vector< Manager >::const_iterator j = v.end();

  i = std::lower_bound( i , j , t , less_cset() );
  j = std::upper_bound( i , j , t , less_cset() );

  std::pair<size_t,size_t> result ;

  result.first  = i - v.begin();
  result.second = j - v.begin();

  return result ;
}

}

//----------------------------------------------------------------------

std::pair<const void*const*,const void*const*>
CSet::p_get( const std::type_info & t ) const
{
  const void * const * b = NULL ;
  const void * const * e = NULL ;

  const std::pair<size_t,size_t> s = span( m_manager , t );

  if ( s.first < s.second ) {
    b = & m_value[ s.first ];
    e = & m_value[ s.second ];
  }
  return SpanVoid( b , e );
}

std::pair<const void*const*,const void*const*>
CSet::p_insert( const Manager & m , const void * v )
{
  std::pair<size_t,size_t> s = span( m_manager , * m.first );
  size_t i ;

  for ( i = s.first ; i < s.second && v != m_value[i] ; ++i );

  if ( i == s.second ) {
    std::vector<Manager>    ::iterator im = m_manager.begin();
    std::vector<const void*>::iterator iv = m_value  .begin();
    std::advance( im , s.second );
    std::advance( iv , s.second );

    m_manager.insert( im , m );
    m_value  .insert( iv , v );
    ++s.second ;
  }
  const void * const * const b = & m_value[ s.first ];
  const void * const * const e = & m_value[ s.second ];
  return SpanVoid( b , e );
}

bool CSet::p_remove( const std::type_info & t , const void * v )
{
  bool result ;

  std::pair<size_t,size_t> s = span( m_manager , t );

  for ( ; s.first < s.second && v != m_value[s.first] ; ++s.first );

  if ( ( result = s.first < s.second ) ) {
    std::vector<Manager>    ::iterator im = m_manager.begin();
    std::vector<const void*>::iterator iv = m_value  .begin();
    std::advance( im , s.first );
    std::advance( iv , s.first );
    m_manager.erase( im );
    m_value  .erase( iv );
  }

  return result ;
}

//----------------------------------------------------------------------

CSet::~CSet()
{
  const size_t n = m_manager.size();
  for ( size_t i = 0 ; i < n ; ++i ) {
    if ( m_manager[i].second ) {
      (*m_manager[i].second)( const_cast<void*>( m_value[i] ) );
    }
  }
}

CSet::CSet() : m_manager(), m_value() {}

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------


#ifdef UNIT_TEST

namespace phdmesh {
namespace unit_test {

class A {
public:
  virtual const char * name() const = 0 ;
  virtual ~A();
};

class B {
public:
  virtual const char * name() const = 0 ;
  virtual ~B();
};

void DoNotDelete( A * a )
{ std::cout << "DoNotDelete(" << a->name() << ")" << std::endl ; }

void DoDelete( A * a )
{
  std::cout << "DoDelete(" << a->name() << ")" << std::endl ;
  delete a ;
}

void DoNotDelete( B * b )
{ std::cout << "DoNotDelete(" << b->name() << ")" << std::endl ; }

void DoDelete( B * b )
{
  std::cout << "DoDelete(" << b->name() << ")" << std::endl ;
  delete b ;
}

class U : public A {
public:
  const char * name() const ;
  ~U() {}
};

class V : public B {
public:
  const char * name() const ;
  ~V() {}
};

class W : public B {
public:
  const char * name() const ;
  ~W() {}
};

class X : public A , public B {
public:
  const char * name() const ;
  ~X() {}
};

class Y : public A , public B {
public:
  const char * name() const ;
  ~Y() {}
};

class Z {
public:
  const char * name() const ;
  ~Z() {}
};

//----------------------------------------------------------------------

int cset()
{
  CSet s ;
  CSet::Span<A> sa ;
  CSet::Span<B> sb ;

  U * u = new U();
  V * v = new V();
  W * w = new W();
  X * x = new X();
  Y * y = new Y();

  std::cout << "s.insert<A>(u,true) " << std::endl ; s.insert<A>(u,true);
  std::cout << "s.insert<B>(v,true) " << std::endl ; s.insert<B>(v,true);
  std::cout << "s.insert<B>(w,true) " << std::endl ; s.insert<B>(w,true);
  std::cout << "s.insert<A>(x) "      << std::endl ; s.insert<A>(x);
  std::cout << "s.insert<B>(x) "      << std::endl ; s.insert<B>(x);

  sa = s.get<A>();
  sb = s.get<B>();

  std::cout
    << "s.get<A>().size() == " << sa.size() << std::endl 
    << "s.get<B>().size() == " << sb.size() << std::endl ;

  for ( unsigned i = 0 ; i < sa.size() ; ++i ) {
    std::cout << "s.get<A>()[" << i << "].name() == "
              << sa[i].name() << std::endl ;
  }
  for ( unsigned i = 0 ; i < sb.size() ; ++i ) {
    std::cout << "s.get<B>()[" << i << "].name() == "
              << sb[i].name() << std::endl ;
  }

  std::cout << "s.remove( <last A> )" << std::endl ; s.remove(& sa.back());

  sa = s.get<A>();
  sb = s.get<B>();

  std::cout << "s.remove( <last B> )" << std::endl ; s.remove(& sb.back());

  sa = s.get<A>();
  sb = s.get<B>();

  delete x ; x = NULL ;

  for ( unsigned i = 0 ; i < sa.size() ; ++i ) {
    std::cout << "s.get<A>()[" << i << "].name() == "
              << sa[i].name() << std::endl ;
  }
  for ( unsigned i = 0 ; i < sb.size() ; ++i ) {
    std::cout << "s.get<B>()[" << i << "].name() == "
              << sb[i].name() << std::endl ;
  }

  std::cout << "s.insert<A>(y,true) " << std::endl ; s.insert<A>(y,true);
  std::cout << "s.insert<B>(y) "      << std::endl ; s.insert<B>(y);

  sa = s.get<A>();
  sb = s.get<B>();

  for ( unsigned i = 0 ; i < sa.size() ; ++i ) {
    std::cout << "s.get<A>()[" << i << "].name() == "
              << sa[i].name() << std::endl ;
  }
  for ( unsigned i = 0 ; i < sb.size() ; ++i ) {
    std::cout << "sb.get<A>()[" << i << "].name() == "
              << sb[i].name() << std::endl ;
  }

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



