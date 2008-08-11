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
 * @date   November 2007
 */

#ifndef util_Span_hpp
#define util_Span_hpp

#include <iterator>

namespace phdmesh {

/** @class Span
 *  Iterate a span of a container defined by begin and end iterators.
 *  Provides forward iterator and const container-like functionality.
 */

template< class IterType ,
          class IterCategory =
             typename std::iterator_traits< IterType >::iterator_category >
class Span ;

//----------------------------------------------------------------------
// Only defined for random access iterators

template< class IterType >
class Span< IterType , std::random_access_iterator_tag > {
public:
  typedef IterType iterator ;

private:
  typedef Span< iterator , std::random_access_iterator_tag > Self ;
  typedef std::iterator_traits< iterator > Traits ;

  iterator m_end ;
  iterator m_iter ;
public:

  //--------------------------------
  // Forward iterator functionality

  typedef typename Traits::value_type      value_type ;
  typedef typename Traits::difference_type difference_type ;
  typedef typename Traits::pointer         pointer ;
  typedef typename Traits::reference       reference ;

  ~Span() {}

  Span() : m_end() { m_iter = m_end ; }

  Span( const Self & rhs ) : m_end( rhs.m_end ) , m_iter( rhs.m_iter ) {}

  Self & operator = ( const Self & rhs )
    { m_end = rhs.m_end ; m_iter = rhs.m_iter ; return *this ; }

  bool operator == ( const Self & rhs ) const
    { return m_end == rhs.m_end && m_iter == rhs.m_iter ; }

  bool operator != ( const Self & rhs ) const
    { return m_end != rhs.m_end || m_iter != rhs.m_iter ; }

  Self & operator ++ () { ++m_iter ; return *this ; }

  Self operator ++ (int) { Self tmp(*this); ++m_iter ; return tmp ; }

  reference operator * ()  const { return *m_iter ; }
  pointer   operator -> () const { return & *m_iter ; }

  //--------------------------------
  // Container-like functionality for random access iterators.

  reference front() const { return *m_iter ; }
  reference back()  const { return m_end[-1] ; }

  iterator begin() const { return m_iter ; }
  iterator end()   const { return m_end ; }

  template<class Iterator>
  Span( Iterator i , Iterator e ) : m_end(e), m_iter(i) {}

  template<class Container>
  explicit
  Span( const Container & c ) : m_end( c.end() ), m_iter( c.begin() ) {}

  template<class Container>
  explicit
  Span( Container & c ) : m_end( c.end() ), m_iter( c.begin() ) {}

  bool empty () const { return ! ( m_iter < m_end ) ; }

  operator bool () const { return m_iter < m_end ; }

  reference operator [] ( difference_type n ) const { return m_iter[n] ; }

  difference_type size() const { return std::distance( m_iter , m_end ); }
};

} // namespace phdmesh

#endif

