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
 * \brief Set of entities of arbitrary types.
 */
class CSet {
public:

  /** \brief  Get member conforming to the given type <b> T </b>.
   *          Return NULL if there is no member of that type.
   */
  template<class T> const T * get() const ;

  /** \brief  Insert a member of a given type <b> T </b>.
   *          Option to transfer ownership of that member.
   *
   *  If a member of the given type already exists then the
   *  insertion operation fails and the existing member is
   *  returned.  If the insertion succeeds then the inserted
   *  member is returned.  For example:
   *
   *  <PRE>
   *    CSet & container = ... ;
   *    const A * const a = new A();
   *    if ( a == container.insert( a ) ) { ... };
   *  </PRE>
   *
   *  If the delete_on_destruction parameter is true then the delete
   *  function will be applied to the inserted member by the
   *  conainter's destructor.
   */
  template<class T>
  const T * insert( const T * , bool delete_on_destruction = false );

  /** \brief  Remove a member of the given type without deleting it.
   *          The caller assumes responsibility for the removed member.
   *          Return if the remove operation was successful.
   */
  template<class T> bool remove( const T * );

  //--------------------------------

  ~CSet();
  CSet();

private:

  typedef void (*DeleteFunction)(void *);

  typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

  const void * p_get( const std::type_info & ) const ;

  const void * p_insert( const Manager & , const void * );

  bool p_remove( const std::type_info & , const void * );

  std::vector< Manager > m_manager ;
  std::vector< const void * > m_value ;

  CSet( const CSet & );
  CSet & operator = ( const CSet & );
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

// Inlined template methods have casting.

namespace phdmesh {

namespace {
template<class T>
void cset_member_delete( void * v ) { delete reinterpret_cast<T*>( v ); }
}

template<class T>
inline
const T * CSet::get() const
{ return (const T*) p_get( typeid(T) ); }

template<class T>
inline
const T * CSet::insert( const T * arg_value , bool delete_on_destruction )
{
  Manager m ;
  m.first = & typeid(T);
  m.second = NULL ;

  if ( delete_on_destruction ) { m.second = & cset_member_delete<T> ; }

  return (const T *) p_insert( m , arg_value );
}

template<class T>
inline
bool CSet::remove( const T * arg_value )
{ return p_remove( typeid(T) , arg_value ); }

} // namespace phdmesh

#endif /* DOXYGEN_COMPILE */

#endif /* util_CSet_hpp */


