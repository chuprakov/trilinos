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

std::ostream & CSet::print( std::ostream & s , const char * separator ) const
{
  static const char space[] = " " ;
  const MemberSet::const_iterator b = m_members.begin();
  const MemberSet::const_iterator e = m_members.end();
  MemberSet::const_iterator i ;

  if ( ! separator ) { separator = space ; }

  for ( i = b ; i != e ; ) {
    if ( i != b ) { s << separator ; }
    CSetMemberBase & m = **i ;
    unsigned count = 0 ;
    while ( i != e && m.m_typeid == (*i)->m_typeid ) {
      ++i ; ++count ;
    }
    s << "typeid(" ;
    s << m.m_typeid.name();
    s << ")#" ;
    s << count ;
  }
  return s ;
}

//----------------------------------------------------------------------

namespace {

// Comparison for sorted vector

struct less_cset {
  bool operator()( const CSetMemberBase * lhs ,
                   const std::type_info & rhs ) const ;
  bool operator()( const std::type_info & lhs ,
                   const CSetMemberBase * rhs ) const ;
};

// On some systems, namely AIX, std::type_info::before(...)
// has a bug where it returns true instead of false for equality.
// Thus we pay a small price on all systems to specifically
// test for and eliminate equality.

bool less_cset::operator()( const CSetMemberBase * lhs ,
                            const std::type_info   & rhs ) const
{ return lhs->m_typeid.before( rhs ) && lhs->m_typeid != rhs ; }

bool less_cset::operator()( const std::type_info   & lhs ,
                            const CSetMemberBase * rhs ) const
{ return lhs.before( rhs->m_typeid ) && lhs != rhs->m_typeid ; }

}

//----------------------------------------------------------------------

Span< CSet::MemberSet::const_iterator >
CSet::span( const std::type_info & t ) const
{
  MemberSet::const_iterator i = m_members.begin();
  MemberSet::const_iterator j = m_members.end();

  i = std::lower_bound( i , j , t , less_cset() );
  j = std::upper_bound( i , j , t , less_cset() );

  return Span< CSet::MemberSet::const_iterator >( i , j );
}

Span< CSet::MemberSet::iterator >
CSet::span( const std::type_info & t )
{
  MemberSet::iterator i = m_members.begin();
  MemberSet::iterator j = m_members.end();

  i = std::lower_bound( i , j , t , less_cset() );
  j = std::upper_bound( i , j , t , less_cset() );

  return Span< CSet::MemberSet::iterator >( i , j );
}

CSet::MemberSet::iterator
CSet::p_insert( CSet::MemberSet::iterator i , CSetMemberBase * v )
{ return m_members.insert( i , v ); }

bool CSet::p_erase( CSet::MemberSet::const_iterator i )
{
  ptrdiff_t n = std::distance( ((const MemberSet &) m_members).begin() , i );

  const bool valid = 0 <= n && n < (ptrdiff_t) m_members.size();

  if ( valid ) { m_members.erase( m_members.begin() + n ); }

  return valid ;
}

//----------------------------------------------------------------------

CSetMemberBase::~CSetMemberBase() {}

CSet::~CSet()
{
  while ( ! m_members.empty() ) {
    delete m_members.back();
    m_members.pop_back();
  }
}

CSet::CSet() : m_members() {}

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
  Span< CSet::iterator<A> > sa ;
  Span< CSet::iterator<B> > sb ;

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
  s.print( std::cout , "  " ) << std::endl ;

  for ( unsigned i = 0 ; i < sa.size() ; ++i ) {
    std::cout << "s.get<A>()[" << i << "].name() == "
              << sa[i].name() << std::endl ;
  }
  for ( unsigned i = 0 ; i < sb.size() ; ++i ) {
    std::cout << "s.get<B>()[" << i << "].name() == "
              << sb[i].name() << std::endl ;
  }

  std::cout << "s.erase( <last A> )" << std::endl ; s.erase(sa.end()-1);

  sa = s.get<A>();
  sb = s.get<B>();

  std::cout << "s.erase( <last B> )" << std::endl ; s.erase(sb.end()-1);

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



