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
 * @date   October 2007
 */

#ifndef util_NamedValue_hpp
#define util_NamedValue_hpp

#include <typeinfo>
#include <utility>
#include <vector>
#include <string>

#include <util/ValueIOS.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** @class NamedValue
 *  @brief ( name , value ) pair
 * 
 *  The public interface for NameValue is essentially the following:
 *
 *    template<typename ValueType>
 *    struct NamedValue {
 *      const std::string name ;
 *            ValueType   value ;
 *      NamedValue( const char * );
 *      NamedValue( const char * , const ValueType & );
 *    };
 *
 *  The ValueType must be default constructable.
 *  A NamedValue must be given a name at construction and
 *  that name is fixed for the lifetime of the NamedValue.
 *  The declarer of a named value owns the object and is
 *  responsible for its destruction.
 *  It is expected that named values will be used where
 *  a regular value would be used, e.g. as a member of a
 *  user's class.
 *
 *    class MyClass {
 *    private:
 *      NamedValue<double> coefficient ;
 *    public:
 *      MyClass() : coefficient("my_coefficient",0.0) {}
 *    };
 *    // ... and elsewhere ...
 *    MyClass * my = new MyClass();
 *    my->coefficient.value = ... ;
 */
template<class Type> class NamedValue ;

//----------------------------------------------------------------------
/** @class NamedValueSet
 *  @brief Ordered set of ( named value , value iostream operator )
 *         that does not own its members.
 *
 *  It is anticipated that a named value will be referenced by
 *  more than one named values set.  As such none of the named
 *  value sets own the name value members.
 *  For example, an IO subsystem could have a set of named value
 *  reference for IO purposes while several other objects use a
 *  subsets of the named values given to the IO subsystem.
 *
 *  It is the responsibility of the declarer of a named value object
 *  to appropriately clean up, i.e. destroy, that object.
 *
 *  The public interface for NamedValueSet is essentially the following:
 *
 *    class NamedValueSet {
 *    public:
 *      NamedValue<T> * get<T>( const std::string & name ,
 *                              const char sep = 0 ) const ;
 *      NamedValue<T> * insert(  NamedValue<T> & val );
 *      NamedValue<T> * replace( NamedValue<T> & val );
 *      void            remove(  NamedValue<T> & val );
 *    };
 *
 *    get<T>( name , sep )
 *      The member of the given value type 'T' and name is returned.
 *      If no such member exists it returns NULL.
 *      If a member exists but of a different type then it throws.
 *      If the separator 'sep' is non zero then 'name' is separated
 *      into sub-names at each occurance of 'sep'.  These sub-names
 *      form a path for searching nested NamedValueSet.
 *
 *    insert( val )
 *      Request to insert 'val' in the container.
 *      If there already exists a member of the given value type and name
 *      then the insert request fails and the existing member is returned.
 *      If the existing member is of a different value type then throws.
 *
 *    replace( val )
 *      Request to insert 'val' in the container.
 *      If there already exists a member of the given value type and name
 *      then the existing member is replaced and returned.
 *      If the existing member is of a different value type then throws.
 *      Effectively (but more efficiently) performs the following:
 *        replace( val )
 *        {
 *          NamedValue<T> * tmp = get<T>( val.name );
 *          if ( tmp ) { remove( *tmp ); }
 *          return insert( val );
 *        }
 *
 *  Advanced interface methods support the following.
 *    ValueIOS<T> objects matched with the members
 *    to enable stream io of an entire NamedValueSet.
 *    Name policy via comparison and good-name functions.
 *
 *  The default comparison is case-insensitive.
 *  The default good-name conforms to the pattern [A-Za-z][_A-Za-z0-9]*.
 */
class NamedValueSet ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Base class to enable NamedValueSet of heterogeneous value types.

template<>
class NamedValue<void> {
public:
  const std::string      name ;
  const std::type_info & value_type ;

  /** Remove this named value from its named value sets.
   *  Virtual to protect against user-code
   *  deleting allocated values via base class.
   */
  virtual ~NamedValue();

  virtual const void * pointer() const = 0 ;
  virtual       void * pointer() = 0 ;

private:
  /** NamedValueSets for which this named value is a member */
  std::vector<NamedValueSet*> m_sets ;

  friend class NamedValueSet ;

  template<typename T> friend class NamedValue ;

  NamedValue( const std::type_info & t , const char * n );

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<void> & );
  NamedValue<void> & operator = ( const NamedValue<void> & );
};

//----------------------------------------------------------------------

template<typename T, unsigned N>
class NamedValue<T[N]> : public NamedValue<void> {
public:
  typedef T          ScalarType ;
  typedef ScalarType ValueType[ N ] ;

  ValueType value ;

  explicit NamedValue( const char * arg_name )
    : NamedValue<void>( typeid(ValueType) , arg_name ) {}

  ~NamedValue() {}

private:

  const void * pointer() const { return & value ; }
        void * pointer()       { return & value ; }

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};


template<typename T>
class NamedValue : public NamedValue<void> {
public:

  typedef T ValueType ;

  ValueType value ;

  explicit NamedValue( const char * arg_name )
    : NamedValue<void>( typeid(ValueType) , arg_name ) {}

  NamedValue( const char * arg_name , const ValueType & arg_val )
    : NamedValue<void>( typeid(ValueType) , arg_name ), value( arg_val ) {}

  ~NamedValue() {}

private:

  const void * pointer() const { return & value ; }
        void * pointer()       { return & value ; }

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
class ValueIOS<NamedValueSet> : public ValueIOS<void> {
public:
  typedef NamedValueSet ValueType ;
  void tell( std::ostream & , unsigned , const ValueType & ) const ;
  void put(  std::ostream & , unsigned , const ValueType & ) const ;
  void get(  std::istream & , ValueType & ) const ;
  ~ValueIOS();
  ValueIOS();
  static const ValueIOS<ValueType> & singleton();
private:
  void tellp( std::ostream & , unsigned , const void * ) const ;
  void putp(  std::ostream & , unsigned , const void * ) const ;
  void getp(  std::istream & , void * ) const ;
};

inline
std::ostream & operator << ( std::ostream & s , const NamedValueSet & v )
{
  ValueIOS<NamedValueSet>::singleton().put( s , 2 , v );
  return s ;
}

inline
std::istream & operator >> ( std::istream & s , NamedValueSet & v )
{
  ValueIOS<NamedValueSet>::singleton().get( s , v );
  return s ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

class NamedValueSet {
public:

  typedef bool (* GoodFunction )( const char * );
  typedef int  (* CompareFunction )( const char * , const char * );
  typedef std::pair< NamedValue<void> * , const ValueIOS<void> * > MemberVoid ;

  const GoodFunction    m_good ;
  const CompareFunction m_compare ;

private:

  std::vector< MemberVoid > m_values ;

  MemberVoid m_get( const std::type_info & ,
                    const std::string & ,
                    const char ) const ;
  MemberVoid m_insert(  MemberVoid );
  MemberVoid m_replace( MemberVoid );

  NamedValueSet & operator = ( const NamedValueSet & );

public:

  //----------------------------------

  ~NamedValueSet();
  NamedValueSet();
  NamedValueSet( const NamedValueSet & );
  explicit NamedValueSet( CompareFunction , GoodFunction = NULL );

  void clear();

  void assign( const std::vector< MemberVoid > & );

  const std::vector< MemberVoid > & get_all() const { return m_values ; }

  //----------------------------------
  // Interface with a default ValueIOS object.

  template<class T>
  NamedValue<T> * get( const std::string & s , const char sep = 0 ) const
    {
      const MemberVoid tmp( m_get( typeid(T) , s , sep ) );
      return static_cast< NamedValue<T> * >( tmp.first );
    }

  template<class T>
  NamedValue<T> * insert( NamedValue<T> & v )
    {
      typedef typename NamedValue<T>::ValueType ValueType ;
      typedef ValueIOS<ValueType> io_type ;
      const io_type & io = io_type::singleton();
      const MemberVoid tmp( m_insert( MemberVoid( & v , & io ) ) );
      return static_cast< NamedValue<T> * >( tmp.first );
    }

  template<class T>
  NamedValue<T> * replace( NamedValue<T> & v )
    {
      typedef ValueIOS<T> io_type ;
      const io_type & io = io_type::singleton();
      const MemberVoid tmp( m_insert( MemberVoid( & v , & io ) ) );
      return static_cast< NamedValue<T> * >( tmp.first );
    }

  /** Remove the reference to the given member from this set. */
  void remove( NamedValue<void> & );

  //----------------------------------
  // Interface with an explicit ValueIOS object, NULL is acceptable.

  template<class T>
  std::pair<NamedValue<T>*,const ValueIOS<T>*>
  get_ios( const std::string & s , const char sep = 0 ) const
    {
      const MemberVoid tmp( m_get( typeid(T) , s , sep ) );
      return std::pair<NamedValue<T>*,const ValueIOS<T>*>(
               static_cast< NamedValue<T> * >( tmp.first ) ,
               static_cast< const ValueIOS<T> * >( tmp.second ) );
    }

  template<class T>
  std::pair< NamedValue<T> * , const ValueIOS<T> * >
  insert_ios( NamedValue<T> & v , const ValueIOS<T> * io )
    {
      const MemberVoid tmp( m_insert( MemberVoid( & v , io ) ) );
      return std::pair<NamedValue<T>*,const ValueIOS<T>*>(
               static_cast< NamedValue<T> * >( tmp.first ) ,
               static_cast< const ValueIOS<T> * >( tmp.second ) );
    }

  template<class T>
  std::pair<NamedValue<T>*,const ValueIOS<T>*>
  replace_ios( NamedValue<T> & v , const ValueIOS<T> * io )
    {
      const MemberVoid tmp( m_replace( MemberVoid( & v , io ) ) );
      return std::pair<NamedValue<T>*,const ValueIOS<T>*>(
               static_cast< NamedValue<T> * >( tmp.first ) ,
               static_cast< const ValueIOS<T> * >( tmp.second ) );
    }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace phdmesh

#endif


