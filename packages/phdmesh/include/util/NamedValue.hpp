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
#include <limits>
#include <utility>
#include <iosfwd>
#include <vector>
#include <string>

#include <util/ValueIOS.hpp>

namespace phdmesh {

class NamedValueSet ;

std::ostream & operator << ( std::ostream & s , const NamedValueSet & v );
std::istream & operator >> ( std::istream & s , NamedValueSet & v );

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
template<class Type = void> class NamedValue ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** @class NamedValueSet
 *  @brief Ordered set of named values that does not own its members.
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
 */
class NamedValueSet {
public:

  typedef int  (* CompareFunction )( const char * , const char * );
  typedef bool (* GoodFunction )( const char * );

  /* Policy for comparing names, default is case-insensitive. */
  const CompareFunction m_compare ;

  /* Policy for a good name, default is [A-Za-z][_A-Za-z0-9]*. */
  const GoodFunction m_good ;

  /* Policy for scalar strean io */
  const ValueIOSPolicy & m_ios ;

  //----------------------------------
  /** The member of the given scalar type and name is returned.
   *  If no such member exists it returns NULL.
   *  If a member exists but of a different type then it throws.
   *  If the separator 'sep' is non zero then 'name' is separated
   *  into sub-names at each occurance of 'sep'.  These sub-names
   *  form a path for searching nested NamedValueSet.
   */
  template<class ScalarType>
  NamedValue<> * find( const std::string & s , const char sep = 0 ) const
    { return m_find( typeid(ScalarType) , s , sep ); }

  //----------------------------------
  /** Request to insert a named value in the container.
   *  If there already exists a member of the given scalar type and name
   *  then the insert request fails and the existing member is returned.
   *  If the existing member is of a different scalar type then throws.
   */
  NamedValue<> * insert( NamedValue<> & v )
    { return m_insert( v , false ); }

  //----------------------------------
  /** Request to insert a named value in the container.
   *  If there already exists a member of the given scalar type and name
   *  then the existing member is replaced and returned.
   *  If the existing member is of a different scalar type then throws.
   */
  NamedValue<> * replace( NamedValue<> & v )
    { return m_insert( v, true ); }

  //----------------------------------
  /** Remove the given member */
  void remove( NamedValue<> & );

  //----------------------------------

  const std::vector<NamedValue<>*> & get_all() const { return m_values ; }

  ~NamedValueSet();
  NamedValueSet();
  NamedValueSet( const NamedValueSet & );

  /** Construct an empty container, with policies.
   *  If any arguments are NULL then the defaults are used.
   */
  NamedValueSet( CompareFunction , GoodFunction , const ValueIOSPolicy * );

  /** Clear the container */
  void clear();

  /** Assign members.  The old members are discarded. */
  void assign( const std::vector<NamedValue<>*> & );

private:

  std::vector< NamedValue<> * > m_values ;

  NamedValue<> * m_find( const std::type_info & ,
                         const std::string & ,
                         const char ) const ;

  NamedValue<> * m_insert( NamedValue<> & , bool );

  NamedValueSet & operator = ( const NamedValueSet & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Base class to enable NamedValueSet of heterogeneous value types.

template<>
class NamedValue<void> {
public:
  const std::string name ;

  template<typename ScalarType>
  void assert_scalar_type() const
    {
      const std::type_info & t = typeid(ScalarType);
      if ( t != scalar_type() ) { throw_type( t ); }
    }

  template<typename ScalarType>
  const ScalarType & get( size_t i = 0 ) const
    {
      assert_scalar_type<ScalarType>();
      return * reinterpret_cast<const ScalarType*>( get_void(i) );
   }

  template<typename ScalarType>
  ScalarType & put( size_t i = 0 )
    {
      assert_scalar_type<ScalarType>();
      return * reinterpret_cast<ScalarType*>( put_void(i) );
    }

  //----------------------------------

  virtual const std::type_info & scalar_type() const = 0 ;

  virtual       size_t get_max() const = 0 ;
  virtual       size_t put_max() const = 0 ;
  virtual       void * put_void( size_t ) = 0 ;
  virtual const void * get_void( size_t ) const = 0 ;

  //----------------------------------
  /** Remove this named value from its named value sets.
   *  Virtual to protect against user-code
   *  deleting allocated values via base class.
   */
  virtual ~NamedValue();

  const std::vector<NamedValueSet*> sets() const { return m_sets ; }

private:
  /** NamedValueSets for which this named value is a member */
  std::vector<NamedValueSet*> m_sets ;

  friend class NamedValueSet ;

  template<typename T> friend class NamedValue ;

  explicit NamedValue( const char * n );

  void throw_type( const std::type_info & ) const ;

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<void> & );
  NamedValue<void> & operator = ( const NamedValue<void> & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue<T&> : public NamedValue<void> {
public:
  typedef T            ScalarType ;
  typedef ScalarType & ValueType ;

  ValueType value ;

  NamedValue( const char * arg_name , ValueType arg_ref )
    : NamedValue<void>(arg_name), value(arg_ref) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return 1 ; }
        size_t put_max() const { return 1 ; }
        void * put_void( size_t i )       { return i ? NULL : & value ; }
  const void * get_void( size_t i ) const { return i ? NULL : & value ; }

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

template<typename T>
class NamedValue<const T &> : public NamedValue<void> {
public:
  typedef T                  ScalarType ;
  typedef const ScalarType & ValueType ;

  ValueType value ;

  NamedValue( const char * arg_name , ValueType arg_ref )
    : NamedValue<void>(arg_name), value(arg_ref) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return 1 ; }
        size_t put_max() const { return 0 ; }
        void * put_void( size_t i )       { return NULL ; }
  const void * get_void( size_t i ) const { return i ? NULL : & value ; }

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue<T*> : public NamedValue<void> {
public:
  typedef T            ScalarType ;
  typedef ScalarType * ValueType ;

  ValueType value ;

  NamedValue( const char * arg_name , ValueType arg_ptr , size_t arg_size )
    : NamedValue<void>(arg_name), value(arg_ptr), m_size( arg_size ) {}

  ~NamedValue() {}

private:

  const size_t m_size ;

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return m_size ; }
        size_t put_max() const { return m_size ; }
        void * put_void(size_t i)       {return i < m_size ? value+i : NULL;}
  const void * get_void(size_t i) const {return i < m_size ? value+i : NULL;}

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue<const T*> : public NamedValue<void> {
public:
  typedef T                  ScalarType ;
  typedef const ScalarType * ValueType ;

  ValueType value ;

  NamedValue( const char * arg_name , ValueType arg_ptr , size_t arg_size )
    : NamedValue<void>(arg_name), value(arg_ptr), m_size( arg_size ) {}

  ~NamedValue() {}

private:

  const size_t m_size ;

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return m_size ; }
        size_t put_max() const { return 0 ; }
        void * put_void(size_t )        {return NULL ; }
  const void * get_void(size_t i) const {return i < m_size ? value+i : NULL;}

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

//----------------------------------------------------------------------

template<typename T, unsigned N>
class NamedValue<T[N]> : public NamedValue<void> {
public:
  typedef T          ScalarType ;
  typedef ScalarType ValueType[ N ] ;

  ValueType value ;

  explicit NamedValue( const char * arg_name ) : NamedValue<void>(arg_name) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return N ; }
        size_t put_max() const { return N ; }
        void * put_void( size_t i )       { return i < N ? & value[i] : NULL ; }
  const void * get_void( size_t i ) const { return i < N ? & value[i] : NULL ; }

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue< std::vector<T> > : public NamedValue<void> {
public:
  typedef T                       ScalarType ;
  typedef std::vector<ScalarType> ValueType ;

  ValueType value ;

  explicit NamedValue( const char * arg_name )
    : NamedValue<void>(arg_name) {}

  explicit NamedValue( const char * arg_name , const ValueType & arg_v )
    : NamedValue<void>(arg_name), value( arg_v ) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

  size_t get_max() const { return value.size(); }
  size_t put_max() const { return std::numeric_limits<size_t>::max(); }

  void * put_void( size_t i )
    { if ( value.size() <= i ) { value.resize(i+1); } return & value[i] ; }

  const void * get_void( size_t i ) const
    { return i < value.size() ? & value[i] : NULL ; }
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue< const std::vector<T> > : public NamedValue<void> {
public:
  typedef T                             ScalarType ;
  typedef const std::vector<ScalarType> ValueType ;

  ValueType value ;

  explicit NamedValue( const char * arg_name , ValueType & arg_v )
    : NamedValue<void>(arg_name), value( arg_v ) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

  size_t get_max() const { return value.size(); }
  size_t put_max() const { return 0 ; }

  void * put_void( size_t ) { return NULL ; }

  const void * get_void( size_t i ) const
    { return i < value.size() ? & value[i] : NULL ; }
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue<const T> : public NamedValue<void> {
public:

  typedef T       ScalarType ;
  typedef const T ValueType ;

  ValueType value ;

  NamedValue( const char * arg_name , const ValueType & arg_val )
    : NamedValue<void>( arg_name ), value( arg_val ) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return 1 ; }
        size_t put_max() const { return 0 ; }
        void * put_void( size_t   )       { return NULL ; }
  const void * get_void( size_t i ) const { return i ? NULL : & value ; }

  // Not implemented:
  NamedValue();
  NamedValue( const NamedValue<ValueType> & );
  NamedValue<ValueType> & operator = ( const NamedValue<ValueType> & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue : public NamedValue<void> {
public:

  typedef T ScalarType ;
  typedef T ValueType ;

  ValueType value ;

  explicit NamedValue( const char * arg_name )
    : NamedValue<void>( arg_name ) {}

  NamedValue( const char * arg_name , const ValueType & arg_val )
    : NamedValue<void>( arg_name ), value( arg_val ) {}

  ~NamedValue() {}

private:

  const std::type_info & scalar_type() const { return typeid(ScalarType); }

        size_t get_max() const { return 1 ; }
        size_t put_max() const { return 1 ; }
        void * put_void( size_t i )       { return i ? NULL : & value ; }
  const void * get_void( size_t i ) const { return i ? NULL : & value ; }

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
  const std::type_info & type() const ;

  ~ValueIOS();
  ValueIOS();

  static const ValueIOSPolicy      & default_policy();
  static const ValueIOS<ValueType> & singleton();
private:
  void tellp( std::ostream & , unsigned , const void * ) const ;
  void putp(  std::ostream & , unsigned , const void * ) const ;
  void getp(  std::istream & , void * ) const ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace phdmesh

#endif


