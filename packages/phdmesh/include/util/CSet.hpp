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
 *  CSet::Span<A> sa = cset.get<A>();
 *  CSet::Span<B> sb = cset.get<B>();
 *
 *  // Remove a member:
 *
 *  {
 *    B * b = ... ;
 *    cset.remove( b ); // Remove never deletes
 *    delete b ;
 *  }
 */
//----------------------------------------------------------------------
/**
 * @class CSet
 * @brief Set of pointers to objects with varied types ordered by type.
 */
class CSet {
public:

  /** @class  Span
   *  @brief  Reference a subset of CSet members of the same type.
   */
  template<class T>
  class Span {
  private:
    friend class CSet ;

    typedef const T * const * iter ;

    iter m_beg ;
    iter m_end ;

    Span( iter b , iter e ) : m_beg(b) , m_end(e) {}

    void validate() { if ( m_end <= m_beg ) { m_end = m_beg = NULL ; } }

  public:

    typedef ptrdiff_t  difference_type ;
    typedef size_t     size_type ;
    typedef const T    value_type ;
    typedef const T *  pointer ;
    typedef const T &  reference ;

    ~Span() {}

    Span() : m_beg(NULL), m_end(NULL) {}

    Span( const Span & rhs ) : m_beg(rhs.m_beg), m_end(rhs.m_end) {}

    Span & operator = ( const Span & rhs )
      { m_beg = rhs.m_beg ; m_end = rhs.m_end ; return *this ; }

    bool operator == ( const Span & rhs ) const
      { return m_beg == rhs.m_beg && m_end == rhs.m_end ; }

    bool operator != ( const Span & rhs ) const
      { return m_beg != rhs.m_beg || m_end != rhs.m_end ; }

    Span & operator ++ () { ++m_beg ; validate(); return *this ; }
    Span & operator -- () { --m_end ; validate(); return *this ; }

    Span operator ++ (int)
      { Span tmp = *this ; ++m_beg ; validate(); return tmp ; }

    Span operator -- (int)
      { Span tmp = *this ; --m_end ; validate(); return tmp ; }

    pointer   operator -> () const { return  *m_beg ; }
    reference operator *  () const { return **m_beg ; }

    reference operator[]( difference_type n ) const
      {
        iter i = n < 0 ? m_end - n : m_beg + n ;
        if ( i < m_beg || m_end <= i ) { i = NULL ; }
        return **i ;
      }

    reference front() const { return **m_beg ; }
    reference back()  const { return *(m_end[-1]); }

    bool empty() const { return ! ( m_beg < m_end ); }

    operator bool () const { return m_beg < m_end ; }

    size_type size() const { return m_end - m_beg ; }
  };

  //--------------------------------
  /** Get members conforming to the given type */
  template<class T> Span<T> get() const ;

  /** Insert and optionally request deletion upon destruction */
  template<class T> Span<T> insert( const T * , bool = false );

  /** Erase a member without deleting */
  template<class T> bool remove( const T * );

  //--------------------------------

  ~CSet();
  CSet();

private:

  typedef void (*DeleteFunction)(void *);

  typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

  typedef std::pair<const void*const*, const void*const*> SpanVoid ;

  SpanVoid p_get( const std::type_info & ) const ;

  SpanVoid p_insert( const Manager & , const void * );

  bool p_remove( const std::type_info & , const void * );

  std::vector< Manager > m_manager ;
  std::vector< const void * > m_value ;

  CSet( const CSet & );
  CSet & operator = ( const CSet & );
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template methods have casting.

namespace phdmesh {

namespace {
template<class T>
void cset_member_delete( void * v ) { delete reinterpret_cast<T*>( v ); }
}

template<class T>
inline
CSet::Span<T> CSet::get() const
{
  const SpanVoid s = p_get( typeid(T) );

  return CSet::Span<T>( reinterpret_cast<const T*const*>( s.first ) ,
                        reinterpret_cast<const T*const*>( s.second ) );
}

template<class T>
inline
CSet::Span<T> CSet::insert( const T * arg_value , bool arg_delete )
{
  Manager m ; m.first = & typeid(T); m.second = NULL ;

  if ( arg_delete ) { m.second = & cset_member_delete<T> ; }

  const SpanVoid s = p_insert( m , arg_value );

  return CSet::Span<T>( reinterpret_cast<const T*const*>( s.first ) ,
                        reinterpret_cast<const T*const*>( s.second ) );
}

template<class T>
inline
bool CSet::remove( const T * arg_value )
{ return p_remove( typeid(T) , arg_value ); }

} // namespace phdmesh

#endif


