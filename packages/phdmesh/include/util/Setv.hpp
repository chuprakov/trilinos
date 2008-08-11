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

#ifndef util_Setv_h
#define util_Setv_h

/*--------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @brief  Enhanced alternative to 'std::set' or 'std::map'.
 *
 *   The 'Setv' alternative to the 'std::set' or 'std::map' template
 *   classes allows explicit memory management of the members of the
 *   container and provides a significantly reduced code-footprint
 *   for the implementation.  Explicit memory management may be
 *   advantageous in complex data structures that have members with
 *   large memory requirements or pointers to members.
 */

// ---------------------------------------------------------------------
// Acknowledgements:
//
//   Most all of the algorithms in this class were obtained from
// the Hewlett-Packard source for the Standard Template Library,
// thus the inclusion of Hewlett-Packard's copyright notice.
// ---------------------------------------------------------------------
/*
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */
/*
Red-black tree class, designed for use in implementing STL
associative containers (set, multiset, map, and multimap).
The insertion and deletion algorithms are based on those in Cormen,
Leiserson, and Rivest, Introduction to Algorithms (MIT Press, 1990),
except that

(1) the header cell is maintained with links not only to the root
but also to the leftmost node of the tree, to enable constant time
begin(), and to the rightmost node of the tree, to enable linear time
performance when used with the generic set algorithms (set_union,
etc.);

(2) when a node being deleted has two children its successor node is
relinked into its place, rather than copied, so that the only
iterators invalidated are those referring to the deleted node.
*/
// ---------------------------------------------------------------------

#include <utility>
#include <iterator>
#include <functional>

namespace phdmesh {

template < typename KeyType > class SetvMember ;

template < class ValueType , bool Forward > class SetvIter ;

template < class ValueType, class KeyCompare, class Allocator > class Setv ;

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template<>
class SetvMember<void> {
public:
  ~SetvMember(); // Destructor removes member from its container
  SetvMember();

  enum { black = 0 , red = 1 };

  SetvMember<void> * parent ; // Binary tree node
  SetvMember<void> * left ;   // Binary tree node
  SetvMember<void> * right ;  // Binary tree node
  unsigned           color ;  // Red-black color

private:
  SetvMember<void>( const SetvMember<void> & );
  SetvMember<void> & operator = ( const SetvMember<void> & );
};

template<bool Forward>
SetvMember<void> * setv_iterate( const SetvMember<void> * );

template<> SetvMember<void> * setv_iterate<true>(  const SetvMember<void> * );
template<> SetvMember<void> * setv_iterate<false>( const SetvMember<void> * );

// ---------------------------------------------------------------------
/** Base class for Setv members.
 *  Objects stored in a Setv container must be derived from this
 *  template base class with the key type for the container.
 *
 *  class MyValueClass : public SetvMember<MyKeyType> { ... };
 *
 *  The key can only be changed by the 'Setv' class.
 */

template < typename KeyType >
class SetvMember : private SetvMember<void> {
public:
  /** Type of key, accessed by Setv container */
  typedef KeyType key_type ;

  /** Query key */
  const key_type & key() const { return m_key ; }

  /** A destroyed node is automatically removed from its container */
  ~SetvMember() {}

  /** Default contructor */
  SetvMember() : SetvMember<void>() {}

  /** Copy constructor */
  SetvMember( const SetvMember<key_type> & rhs )
    : SetvMember<void>(), m_key( rhs.m_key ) {}

  /** explicit Key constructor */
  explicit SetvMember( const KeyType & k )
    : SetvMember<void>(), m_key( k ) {}

private:

  SetvMember<key_type> & operator = ( const SetvMember<key_type> & m );

  key_type m_key ;

  template< class T , class C , class A > friend class Setv ;
  template< class T , bool F > friend class SetvIter ;
};

//----------------------------------------------------------------------
/** @class SetvIter
 *  Template class for the Setv iterator types.
 */
template < class Type , bool Forward>
class SetvIter : public std::iterator<std::bidirectional_iterator_tag,Type>
{
private:

  Type * n ;

  SetvIter( Type * x ) : n(x) {}

  template < class T , bool F >            friend class SetvIter ;
  template < class T , class C , class A > friend class Setv ;

public:

  /** [EXTENSION] Conversion to bool: true if a dereferenceable iterator */
  operator bool () const { return n && n->parent ; }

  // Constructors & assignment operator

  SetvIter() : n(NULL) {}

  /** Copy from other similar iterators.
   *  Compilation will fail when assigning non-const to const.
   */
  template<class T,bool F> SetvIter( const SetvIter<T,F> & x ) : n(x.n) {}

  /** Assign from other similar iterators.
   *  Compilation will fail when assigning non-const to const.
   */
  template<class T,bool F>
    SetvIter<Type,Forward> & operator = ( const SetvIter<T,F> & x )
      { n = x.n ; return *this ; }

  template<class T,bool F>
    bool operator == ( const SetvIter<T,F> & y ) const { return n == y.n ; }

  template<class T,bool F>
    bool operator != ( const SetvIter<T,F> & y ) const { return n != y.n ; }

  Type & operator * () const
    { return *(operator bool() ? n : reinterpret_cast<Type*>(NULL) ); }
  Type * operator ->() const
    { return  (operator bool() ? n : reinterpret_cast<Type*>(NULL) ); }

  SetvIter<Type,Forward> & operator++()
    { n = static_cast<Type*>( setv_iterate<Forward>(n) ); return *this ; }

  SetvIter<Type,Forward> & operator--()
    { n = static_cast<Type*>( setv_iterate<!Forward>(n) ); return *this ; }

  SetvIter<Type,Forward> operator++(int)
    {
      Type * const t = n ;
      n = static_cast<Type*>( setv_iterate<Forward>(n) );
      return SetvIter<Type,Forward>(t);
    }

  SetvIter<Type,Forward> operator--(int)
    {
      Type * const t = n ;
      n = static_cast<Type*>( setv_iterate<!Forward>(n) );
      return SetvIter<Type,Forward>(t);
    }
};

//----------------------------------------------------------------------

template<>
class Setv<void,void,void> : private SetvMember<void> {
protected:
  friend class SetvMember<void> ; // For destructor to remove

  Setv();
  ~Setv();

  size_t size() const { return m_size ; }

  SetvMember<void> * nRoot() const { return m_header.parent ; }

  SetvMember<void> * nEnd() const
    { return const_cast<SetvMember<void>*>( & m_right_end ); }

  SetvMember<void> * nREnd() const
    { return const_cast<SetvMember<void>*>( & m_left_end ); }

  SetvMember<void> * nBegin()  const
    { return ( m_left_end.right != nREnd() ) ? m_left_end.right : nEnd(); }

  SetvMember<void> * nRBegin() const
    { return ( m_right_end.left != nEnd() ) ? m_right_end.left : nREnd(); }

  void initialize();

  void insert( SetvMember<void> * y , SetvMember<void> * z , bool z_lt_y );

  void remove( SetvMember<void> * );

  SetvMember<void> * unbalancing_removal( SetvMember<void> ** n );

  static Setv<void,void,void> * container( const SetvMember<void> * );

private:

  SetvMember<void> & m_header ;    // head node is this
  SetvMember<void>   m_left_end ;  // 'rend()' node
  SetvMember<void>   m_right_end ; // 'end()' node
  size_t             m_size ;

  // root      == m_header.parent
  // leftmost  == m_left_end.right
  // rightmost == m_right_end.left
};

//----------------------------------------------------------------------

template < class ValueType ,
           class KeyCompare = std::less<typename ValueType::key_type> ,
           class Allocator = std::allocator<ValueType> >
class Setv : private Setv<void,void,void> {
public:
  // Types:
  typedef typename ValueType::key_type key_type ;
  typedef ValueType                    value_type ;
  typedef KeyCompare                   key_compare ;
  typedef Allocator                    allocator_type ;

  typedef typename allocator_type::reference       reference ;
  typedef typename allocator_type::const_reference const_reference ;
  typedef typename allocator_type::pointer         pointer ;
  typedef typename allocator_type::const_pointer   const_pointer ;
  typedef typename allocator_type::size_type       size_type ;
  typedef typename allocator_type::difference_type difference_type ;

  // Iterators are bidirectional and have as an extension a
  // conversion-to-bool operator that returns true for a
  // dereferencable iterator.

  typedef SetvIter<      value_type,true>  iterator ;
  typedef SetvIter<const value_type,true>  const_iterator ;
  typedef SetvIter<      value_type,false> reverse_iterator ;
  typedef SetvIter<const value_type,false> const_reverse_iterator ;

  struct value_compare : public std::binary_function<value_type,value_type,bool>
  {
    protected:
      key_compare comp ;
    public:
      bool operator()(const value_type& x, const value_type& y) const {
        return comp( x.SetvMember<key_type>::key() ,
                     y.SetvMember<key_type>::key() );
      }
  };

  // Construct/copy/destroy:

  ~Setv();

  Setv( const key_compare    & arg_compare = key_compare(),
        const allocator_type & arg_alloc = allocator_type() )
    : Setv<void,void,void>(),
      alloc(arg_alloc), key_less(arg_compare), value_less() {}

  Setv( const Setv<value_type,key_compare,allocator_type> & );

  Setv<value_type,key_compare,allocator_type> & operator =
    ( const Setv<value_type,key_compare,allocator_type> & );

  allocator_type get_allocator() const { return alloc ; }

  // Iterators:
  iterator begin() const
    { return iterator(
        static_cast<value_type*>( Setv<void,void,void>::nBegin() ) ); }

  iterator end() const
    { return iterator(
        static_cast<value_type*>( Setv<void,void,void>::nEnd() ) ); }

  reverse_iterator rbegin() const
    { return reverse_iterator(
        static_cast<value_type*>( Setv<void,void,void>::nRBegin() ) ); }

  reverse_iterator rend() const
    { return reverse_iterator(
        static_cast<value_type*>( Setv<void,void,void>::nREnd() ) ); }

  // Capacity:
  bool      empty() const { return Setv<void,void,void>::size() == 0 ; }
  size_type size()  const { return Setv<void,void,void>::size() ; }
  size_type max_size() const { return alloc.max_size(); }

  // Modifiers (single member only):

  /** [EXTENSION] Inserts the given value with key() == key
   *
   *  If return.second == true then the value was removed
   *  from its existing container, had its key() set to 'key',
   *  and was inserted into this container.  If value == NULL
   *  then a value is allocated.
   *
   *  If return.second == false then an existing value
   *  was found and returned, and the input value is untouched.
   *
   *  The value must be individually deallocatable via the allocator_type,
   *  i.e. allocator.deallocate(value,1) is valid.
   */
  std::pair<iterator,bool>
    insert( const key_type & key , value_type * value = NULL );

  /** [EXTENSION]  Inserts the given value with its key already set */
  std::pair<iterator,bool> insert( value_type * );

  /** Retrieve or create a value with the given key */
  value_type & operator[]( const key_type & key );

  /** Destroys member */
  void erase( iterator position );

  /** [EXTENSION] Destroys member */
  void erase( value_type * );

  /** Destroys member */
  size_type erase( const key_type & );

  /** [EXTENSION] Removes member from container but does not deallocate
   *  The caller assumes responsibility for the value.
   */
  void remove( value_type & v )
    { Setv<void,void,void>::remove( &v ); }

  // Modifiers: whole container:

  // void swap( Setv<value_type,key_compare,allocator_type> & );

  /** Destroys all member according to allocator_type */
  void clear();

  // Observers:
  key_compare   key_comp()   const { return key_less ; }
  value_compare value_comp() const { return value_less ; }

  // Set operations:
  iterator  find(  const key_type & ) const ;
  size_type count( const key_type & ) const ;

  iterator  lower_bound( const key_type & ) const ;
  iterator  upper_bound( const key_type & ) const ;

  /** [EXTENSION] Verify that iteration of the container is properly ordered */
  bool verify_ordering() const ;

  /** [EXTENSION] Return the container of this member [extension] */
  static
  Setv<value_type,key_compare,allocator_type> * container( const value_type & );

private:

  typedef SetvMember<key_type> MemberType ;

  allocator_type alloc ;
  key_compare    key_less ;
  value_compare  value_less ;
};

}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Implementation from here forward

namespace phdmesh {

//----------------------------------------------------------------------

template<class T,class C,class M>
inline
Setv<T,C,M>::~Setv()
{ clear(); }

template<class T,class C,class M>
inline
std::pair<typename Setv<T,C,M>::iterator,bool>
Setv<T,C,M>::insert( typename Setv<T,C,M>::value_type * v )
{
  bool flag = false ;

  if ( v != NULL ) {
    const typename Setv<T,C,M>::key_type & k = v->key();

    phdmesh::SetvMember<void> * y = nEnd();
    phdmesh::SetvMember<void> * x = nRoot();

    flag = true ;

    while ( x ) {
      y = x ;
      x = ( flag = key_less( k , static_cast<MemberType*>(x)->key() ) )
        ? x->left : x->right ;
    }

    /* flag = k < y , check previous value if exists */

    const bool k_lt_y = flag ;

    x = flag && y != nBegin() ? ( flag = false , setv_iterate<false>(y) ) : y ;

    if ( flag || ( flag = key_less(static_cast<MemberType*>(x)->key(),k) ) ) {
      x = v ;
      Setv<void,void,void>::insert( y , x , k_lt_y );
    }
    else {
      v = static_cast<typename Setv<T,C,M>::value_type*>( x );
    }
  }
  return std::pair<typename Setv<T,C,M>::iterator,bool>( iterator(v), flag );
}

template<class T,class C,class M>
inline
std::pair<typename Setv<T,C,M>::iterator,bool>
Setv<T,C,M>::insert( const typename Setv<T,C,M>::key_type   & k ,
                           typename Setv<T,C,M>::value_type * v )
{
  phdmesh::SetvMember<void> * y = nEnd();
  phdmesh::SetvMember<void> * x = nRoot();

  bool flag = true ;

  while ( x ) {
    y = x ;
    x = ( flag = key_less( k , static_cast<MemberType*>(x)->key() ) )
      ? x->left : x->right ;
  }

  /* flag = k < y , check previous value if exists */

  const bool k_lt_y = flag ;

  x = flag && y != nBegin() ? ( flag = false , setv_iterate<false>(y) ) : y ;

  if ( flag || ( flag = key_less(static_cast<MemberType*>(x)->key(),k) ) ) {
    if ( v == NULL ) {
      v = alloc.allocate(1);
      new(v) value_type();
    }
    v->SetvMember<key_type>::m_key = k ;
    x = v ;
    Setv<void,void,void>::insert( y , x , k_lt_y );
  }
  else {
    v = static_cast<typename Setv<T,C,M>::value_type*>( x );
  }
  return std::pair<typename Setv<T,C,M>::iterator,bool>( iterator(v), flag );
}

template<class T,class C,class M>
inline
typename Setv<T,C,M>::value_type &
Setv<T,C,M>::operator[]( const key_type & k )
{
  std::pair<typename Setv<T,C,M>::iterator,bool> result = insert(k);
  return *result.second ;
}

template<class T,class C,class M>
inline
void Setv<T,C,M>::erase( typename Setv<T,C,M>::value_type * p )
{
  Setv<void,void,void>::remove( p );
  alloc.destroy(    p );
  alloc.deallocate( p , 1 );
}

template<class T,class C,class M>
inline
void Setv<T,C,M>::erase( typename Setv<T,C,M>::iterator p )
{
  if ( p.operator bool() ) { erase( p.n ); }
}

template<class T,class C,class M>
inline
typename Setv<T,C,M>::size_type
Setv<T,C,M>::erase( const typename Setv<T,C,M>::key_type & k )
{
  iterator i = find(k);
  return i != end() ? ( erase( i ) , 1 ) : 0 ;
}

template<class T,class C,class M>
void Setv<T,C,M>::clear()
{
  if ( Setv<void,void,void>::size() ) {
    phdmesh::SetvMember<void> * n = nBegin();
    phdmesh::SetvMember<void> * t ;
    while ( ( t = Setv<void,void,void>::unbalancing_removal( &n ) ) ) {
      alloc.destroy(    static_cast<value_type*>(t) );
      alloc.deallocate( static_cast<value_type*>(t) , 1 );
    }
  }
  Setv<void,void,void>::initialize();
}

template<class T,class C,class M>
inline
typename Setv<T,C,M>::iterator  
Setv<T,C,M>::lower_bound( const typename Setv<T,C,M>::key_type & k ) const  
{
  phdmesh::SetvMember<void> * y = nEnd();
  phdmesh::SetvMember<void> * x = nRoot();
  while ( x ) x = key_less(static_cast<MemberType*>(x)->key(),k)
                ? x->right : ( y = x )->left ;
  return iterator( static_cast<value_type*>(y) );
}

template<class T,class C,class M>
inline
typename Setv<T,C,M>::iterator  
Setv<T,C,M>::upper_bound( const typename Setv<T,C,M>::key_type & k ) const  
{
  phdmesh::SetvMember<void> * y = nEnd();
  phdmesh::SetvMember<void> * x = nRoot();
  while ( x ) x = key_less(k,static_cast<MemberType*>(x)->key())
                ? ( y = x )->left : x->right ;
  return iterator( static_cast<value_type*>(y) );
}

template<class T,class C,class M>
inline
typename Setv<T,C,M>::iterator  
Setv<T,C,M>::find( const typename Setv<T,C,M>::key_type & k ) const  
{
  typename Setv<T,C,M>::iterator i = lower_bound(k); // k <= i->key();

  // If 'end()' or k < y->key() then not found
  if ( i != end() && key_less(k,(*i).MemberType::key()) ) i = end();

  return i ;
}

template<class T,class C,class M>
inline
typename Setv<T,C,M>::size_type
Setv<T,C,M>::count( const typename Setv<T,C,M>::key_type & k ) const  
{
  typename Setv<T,C,M>::iterator i = lower_bound(k); // k <= i->key();
  return i == end() || key_less(k,(*i).SetvMember<T>::key()) ? 0 : 1 ;
}

template<class T,class C,class M>
inline
Setv<T,C,M>::Setv( const Setv<T,C,M> & V )
  : Setv<void,void,void>()
{
  for ( const_iterator i = V.begin() ; i != V.end() ; ++i ) {
    pointer a = alloc.allocate(1);
    alloc.construct(a,*i);
    insert( i->SetvMember<key_type>::key() , a );
  }
}

template<class T,class C,class M>
inline
Setv<T,C,M> &
Setv<T,C,M>::operator = ( const Setv<T,C,M> & V )
{
  clear();
  for ( const_iterator i = V.begin() ; i != V.end() ; ++i ) {
    pointer a = alloc.allocate(1);
    alloc.construct(a,*i);
    insert( i->SetvMember<key_type>::key() , a );
  }
}

  /** [EXTENSION] Return the container of this member [extension] */
template<class T,class C,class M>
inline
Setv<T,C,M> *
Setv<T,C,M>::container( const typename Setv<T,C,M>::value_type & v )
{
  Setv<void,void,void> * const c = Setv<void,void,void>::container(&v);
  return c ? static_cast<Setv<T,C,M>*>( c ) :
             reinterpret_cast<Setv<T,C,M>*>( NULL );
}

template<class T,class C,class M>
inline
bool Setv<T,C,M>::verify_ordering() const
{
  size_type count = 0 ;
  iterator  i = begin();
  iterator  j ;

  while ( i != end() &&
         ( ++(j=i) == end() || key_less( i->MemberType::key() ,
                                         j->MemberType::key() ) )
         && --j == i )
    { ++i ; ++count ; }

  return ( i == end() && count == size() );
}

// ---------------------------------------------------------------------

}

#endif

