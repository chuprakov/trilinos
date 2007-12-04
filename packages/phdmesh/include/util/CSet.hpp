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
 * @date   November 2006
 */

#ifndef util_CSet_hpp
#define util_CSet_hpp

#include <typeinfo>
#include <vector>
#include <iosfwd>

#include <util/Span.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/**
 * @class CSet
 * @brief Multiset of entities of arbitrary types.
 *
 *  Example usage of the three methods:
 *
 *  class A { ... };
 *  class B { ... };
 *
 *  CSet cset ;
 *
 *  // Insert pointers to objects:
 *
 *  cset.insert<A>( new A ); // Do not delete on destruction
 *  cset.insert<B>( new B , true ); // Delete on destruction
 *  cset.insert<B>( new B , true ); // Delete on destruction
 *
 *  // Query the collection of objects of a given type:
 *
 *  Span< CSet::iterator<A> > sa = cset.get<A>();
 *  Span< CSet::iterator<B> > sb = cset.get<B>();
 *
 *  // Erase a member:
 *
 *  {
 *    B * b = cset.erase( sb.begin() ); // Erase never deletes
 *    delete b ;
 *  }
 */
//----------------------------------------------------------------------

class CSet ;

class CSetMemberBase {
public:
  virtual ~CSetMemberBase();
  const std::type_info & m_typeid ;
protected:
  bool m_delete ;
  explicit CSetMemberBase( const std::type_info & t )
    : m_typeid(t), m_delete(false) {}
private:
  CSetMemberBase();
  CSetMemberBase( const CSetMemberBase & );
  CSetMemberBase & operator = ( const CSetMemberBase & );
  friend class CSet ;
};

namespace {

template<class T>
class CSetMember : public CSetMemberBase {
public:
  T * m_value ;
  virtual ~CSetMember();
  explicit CSetMember( T * v ) : CSetMemberBase( typeid(T) ), m_value(v) {}
private:
  CSetMember();
  CSetMember( const CSetMember<T> & );
  CSetMember<T> & operator = ( const CSetMember<T> & );
};

template<class T>
CSetMember<T>::~CSetMember() { if ( m_delete ) { delete m_value ; } }

}

//----------------------------------------------------------------------
/**
 * @class CSet
 * @brief Set of pointers to objects with varied types.
 */
class CSet {
public:

  typedef std::vector< CSetMemberBase * > MemberSet ;

  template<class T>
  class iterator : public std::iterator<std::random_access_iterator_tag,T> {
  private:
    friend class CSet ;

    MemberSet::const_iterator i ;

    T * value( ptrdiff_t n ) const
      { return static_cast<CSetMember<T>*>(*(i+n))->m_value ; }

    T * value() const { return static_cast<CSetMember<T>*>(*i)->m_value ; }

  public:

    iterator( MemberSet::iterator rhs ) : i(rhs) {}
    iterator( MemberSet::const_iterator rhs ) : i(rhs) {}

    ~iterator() {}
    iterator() {}
    iterator( const iterator<T> & rhs ) : i( rhs.i ) {}
    iterator<T> & operator = ( const iterator<T> & rhs )
      { i = rhs.i ; return *this ; }

    T & operator[]( ptrdiff_t n ) const { return *value(n); }
    T & operator *  () const { return *value(); }
    T * operator -> () const { return value(); }

    bool operator == ( const iterator<T> & rhs ) const { return i == rhs.i ; }
    bool operator != ( const iterator<T> & rhs ) const { return i != rhs.i ; }
    bool operator <  ( const iterator<T> & rhs ) const { return i <  rhs.i ; }
    bool operator >  ( const iterator<T> & rhs ) const { return i >  rhs.i ; }
    bool operator <= ( const iterator<T> & rhs ) const { return i <= rhs.i ; }
    bool operator >= ( const iterator<T> & rhs ) const { return i >= rhs.i ; }

    iterator<T> & operator ++ () { ++i ; return *this ; }
    iterator<T> & operator -- () { --i ; return *this ; }
    iterator<T>   operator ++ (int) { return iterator<T>( i++ ); }
    iterator<T>   operator -- (int) { return iterator<T>( i-- ); }
    iterator<T> & operator += (ptrdiff_t n) { i+=n ; return *this; }
    iterator<T> & operator -= (ptrdiff_t n) { i-=n ; return *this; }
    iterator<T>   operator +  (ptrdiff_t n) const {return iterator<T>(i+n);}
    iterator<T>   operator -  (ptrdiff_t n) const {return iterator<T>(i-n);}

    ptrdiff_t operator - ( const iterator<T> & rhs ) const
      { return i - rhs.i ; }
  };

  //--------------------------------
  /** Insert and optionally request deletion upon destruction */
  template<class T> iterator<T> insert( T * , bool = false );

  /** Erase a member without deleting */
  template<class T> T * erase( iterator<T> );

  /** Get members conforming to the given type */
  template<class T> Span< iterator<T> > get() const ;

  //--------------------------------

  ~CSet();
  CSet();

  std::ostream & print( std::ostream & , const char * separator ) const ;

private:

  Span< MemberSet::const_iterator > span( const std::type_info & ) const ;
  Span< MemberSet::iterator >       span( const std::type_info & );

  MemberSet::iterator p_insert( MemberSet::iterator , CSetMemberBase * );

  bool p_erase( MemberSet::const_iterator );

  MemberSet m_members ;
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template methods have casting.

namespace phdmesh {

template<class T>
inline
Span< CSet::iterator<T> > CSet::get() const
{
  const Span< MemberSet::const_iterator > s = span( typeid(T) );
  return Span< CSet::iterator<T> >( s.begin() , s.end() );
}

template<class T>
inline
CSet::iterator<T> CSet::insert( T * arg_value , bool arg_delete )
{
  Span< MemberSet::iterator > s = span( typeid(T) );

  for ( ; s && arg_value != static_cast<CSetMember<T>*>(*s)->m_value ; ++s );

  const MemberSet::iterator i =
    s ? s.begin() : p_insert( s.end() , new CSetMember<T>(arg_value) );

  (*i)->m_delete = arg_delete ;

  return CSet::iterator<T>( i );
}

template<class T>
inline
T * CSet::erase( CSet::iterator<T> j )
{
  T * const p = j.value();
  return p_erase( j.i ) ? p : (T*) NULL ;
}

} // namespace phdmesh

#endif


