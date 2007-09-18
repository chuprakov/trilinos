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
#include <utility>
#include <vector>
#include <string>

namespace phdmesh {

//----------------------------------------------------------------------
/**
 * @class CSet
 * @brief Multiset of entities with varying virtual base classes
 *
 * @par Simple inheritence usage
 *
 *    Let 'A' and 'B' be virtual base classes and
 *    let 'U', 'V', and 'W' be implementation classes.
 *
 *      class A : public phdmesh::CSetMember<A> { ... };
 *      class B : public phdmesh::CSetMember<B> { ... };
 *
 *    The 'phdmesh::CSetMember<class T>' base class enables
 *    implementations of 'A' and 'B' to be inserted into a
 *    CSet container.  DO NOT use virtual inheritence in the
 *    above derivations.
 *
 *      class U : public A { ... };
 *      class V : public B { ... };
 *      class W : public B { ... };
 *
 *    Insert implementations of U, V, and W into a virtual set.
 *
 *      CSet s ;
 *      U * inst_u = new U();
 *      V * inst_v = new V();
 *      W * inst_w = new W();
 *
 *      s.insert<A>( inst_u ) == 0 // The first  member of class 'A'
 *      s.insert<B>( inst_v ) == 0 // The first  member of class 'B'
 *      s.insert<B>( inst_w ) == 1 // The second member of class 'B'
 *
 *      s.count<A>() == 1
 *      s.count<B>() == 2
 *
 *    The default ordinal to the 'get' method is zero
 *
 *      A * au = s.get<A>();   au == inst_u
 *      B * bv = s.get<B>();   bv == inst_v
 *      B * bw = s.get<B>(1);  bw == inst_w
 *
 * @par Geting all implementations of a given base class
 *
 *      std::vector<const A *> a_all = s.get_all<A>();
 *
 * @par Destruction of member implementations
 *
 *    The 'insert' method has a second argument 'destroy' which
 *    is a function to destroy the inserted member when ever the
 *    container is destroyed.  The default value for 'destroy' is
 *    a function that calls the 'delete' operator on the member.
 *    If the 'insert' method is given 'NULL' then nothing is done
 *    to the inserted member when the container is destroyed.
 *
 *      U & inst_u = singleton_of_u();
 *
 *      s.insert<A>( inst_u , NULL ); // Do not call 'delete' on 'inst_u'
 *
 * @par Multiple interitence usage
 *
 *    Let 'A' and 'B' be virtual base classes as before,
 *    let 'U', 'V', and 'W' be implementation classes as before, and
 *    let 'X' implement both 'A' and 'B'
 *
 *      class X : public A , public B { ... }
 *
 *    Insert impelementations of U, V, and W as before.
 *    Insert X into the same virtual set, but only destroy it once!
 *
 *      X * inst_x = new X();
 *
 *      s.insert<A>( inst_x ) == 1        // The second member of class 'A'
 *      s.insert<B>( inst_x , NULL ) == 2 // The third  member of class 'B'
 *
 *      s.count<A>() == 2
 *      s.count<B>() == 3
 *
 *      A * ax = s.get<A>(1);   ax == inst_x
 *      B * bx = s.get<B>(2);   bx == inst_x
 *
 */
//----------------------------------------------------------------------

class CSet ;
template<class Type> struct CSetMember ;

//----------------------------------------------------------------------
/**
 * @class CSet
 * @brief Set of entities with varied virtual base classes.
 */

class CSet {
public:

  static void default_destroy( const CSetMember<void> * m );

  typedef void (*MemberDestroy)( const CSetMember<void> * );

  //--------------------------------

  ~CSet();
  CSet();

  /** Insert a member of the given base type.
   *  Return the ordinal of the member among members with the same base type.
   *  Ordinals are incremented with each insertion of the same base type.
   */

  template<class T>
    unsigned insert( const T * , MemberDestroy = & default_destroy );

  /** Replace a member of the given base type.
   *  Return ordinal of the replaced member.
   *  The caller assumes ownership of the replaced member.
   */

  template<class T>
    unsigned replace( const T * c_old ,
                      const T * c_new , MemberDestroy = & default_destroy );

  /** Get count of members of the given base type */

  template<class T> unsigned count() const ;

  /** Get a member of the given base type with the given ordinal.  */

  template<class T>
    const typename T::CSet_Type * get( unsigned ordinal = 0 ) const ;

  /** Get all members of a given base type.  */

  template<class T>
    std::vector<const typename T::CSet_Type *> get_all() const ;

  void print( std::string & , const char * separator ) const ;

private:

  typedef std::pair< const CSetMember<void> * ,
                     void (*)( const CSetMember<void> * ) > ValueType ;

  typedef std::pair< std::vector<ValueType>::const_iterator ,
                     std::vector<ValueType>::const_iterator > BoundsType ;

  static
  std::vector<ValueType>::const_iterator
  m_upper_bound( std::vector<ValueType>::const_iterator ,
                 std::vector<ValueType>::const_iterator ,
                 const std::type_info & );
  static 
  std::vector<ValueType>::const_iterator
  m_lower_bound( std::vector<ValueType>::const_iterator ,
                 std::vector<ValueType>::const_iterator ,
                 const std::type_info & );

  static
  std::vector<ValueType>::iterator
  m_lower_bound( std::vector<ValueType>::iterator ,
                 std::vector<ValueType>::iterator ,
                 const std::type_info & );

  BoundsType m_bounds( const std::type_info & ) const ;

  unsigned m_insert( const CSetMember<void> * , MemberDestroy ,
                     const char * const );

  unsigned m_replace( const CSetMember<void> * ,
                      const CSetMember<void> * , MemberDestroy ,
                      const char * const );

  std::vector< ValueType > m_inst ;
};

//----------------------------------------------------------------------

template<>
struct CSetMember<void> {
  /** The implementation class may provide a name */
  virtual const char * name() const ;

  /** Virtual destructor for destruction of the CSet container. */
  virtual ~CSetMember() {}

  /** The base class type (container key) for this member */
  const std::type_info & m_cset_type ;

private:

  CSetMember( const std::type_info & t ) : m_cset_type(t) {}

  CSetMember();
  CSetMember( const CSetMember<void> & );
  CSetMember<void> & operator = ( const CSetMember & );

  template<class T> friend class CSetMember ;
};

template<class T>
struct CSetMember : public CSetMember<void> {
  typedef T CSet_Type ;

  virtual ~CSetMember() {}

  CSetMember()                     : CSetMember<void>( typeid(CSet_Type) ) {}
  CSetMember(const CSetMember<T> &) : CSetMember<void>( typeid(CSet_Type) ) {}
  CSetMember<T> & operator = ( const CSetMember<T> & ) {}
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template methods have casting.

namespace phdmesh {

template<class T>
inline
unsigned CSet::count() const
{
  const BoundsType b = m_bounds( typeid(typename T::CSet_Type) );
  return b.second - b.first ;
}

template<class T>
inline
const typename T::CSet_Type * CSet::get( unsigned ordinal ) const
{
  BoundsType b = m_bounds( typeid(typename T::CSet_Type) );

  return ( b.first += ordinal ) < b.second
         ? static_cast<const typename T::CSet_Type *>( b.first->first )
         : reinterpret_cast<const typename T::CSet_Type*>( NULL );
}

template<class T>
inline
std::vector<const typename T::CSet_Type *> CSet::get_all() const
{
  std::vector<const typename T::CSet_Type *> result ;

  BoundsType b = m_bounds( typeid(typename T::CSet_Type) );

  result.reserve( b.second - b.first );

  for ( ; b.first != b.second ; ++b.first ) {
    result.push_back(
      static_cast<const typename T::CSet_Type *>(b.first->first) );
  }
  return result ;
}

template<class T>
inline
unsigned CSet::insert( const T * t , CSet::MemberDestroy destroy )
{
  // Require the proper inheritence at insertion.

  const CSetMember<typename T::CSet_Type> * const ts = t ;
  const CSetMember<void>                  * const tv = ts ;

  return m_insert( tv , destroy , typeid(T).name() );
}


template<class T>
inline
unsigned CSet::replace( const T * c_old ,
                        const T * c_new , CSet::MemberDestroy destroy )
{
  // Require the proper inheritence at removal.

  const CSetMember<typename T::CSet_Type> * const ts_old = c_old ;
  const CSetMember<void>                  * const tv_old = ts_old ;
  const CSetMember<typename T::CSet_Type> * const ts_new = c_new ;
  const CSetMember<void>                  * const tv_new = ts_new ;

  return m_replace( tv_old , tv_new , destroy , typeid(T).name() );
}


} // namespace phdmesh

#endif


