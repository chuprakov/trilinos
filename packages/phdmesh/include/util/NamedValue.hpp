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
#include <iosfwd>
#include <vector>
#include <string>

#include <util/Reference.hpp>
#include <util/ValueIO.hpp>

namespace phdmesh {

class NamedValueSet ;

std::ostream & operator << ( std::ostream & , const NamedValueSet & );
std::istream & operator >> ( std::istream & ,       NamedValueSet & );

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

template<>
struct ValueIO<NamedValueSet> {
  enum { binary_size = 0 };
  static void   write( std::ostream &, const NamedValueSet*, size_t );
  static size_t read(  std::istream &,       NamedValueSet*, size_t );
  static void   tell(  std::ostream &, const NamedValueSet*, size_t, size_t );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Base class to enable NamedValueSet of heterogeneous value types.

template<>
class NamedValue<void> {
public:
  const std::string       name ;
  const Reference<void> & reference ;

  //----------------------------------
  /** Remove this named value from its named value sets.
   *  Virtual to protect against user-code
   *  deleting allocated values via base class.
   */
  virtual ~NamedValue();

  const std::vector<NamedValueSet*> sets() const { return m_sets ; }

protected:

  NamedValue( const char * , const Reference<void> & );

private:
  /** NamedValueSets for which this named value is a member */
  std::vector<NamedValueSet*> m_sets ;

  friend class NamedValueSet ;

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------

template<typename T, unsigned N>
class NamedValue<T[N]> : public NamedValue<void> {
public:

  explicit NamedValue( const char * arg_name )
    : NamedValue<void>(arg_name,m_ref), m_ref( value , N ) {}

  ~NamedValue() {}

  T value[N] ;

private:

  Reference<T*> m_ref ;

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue<T*> : public NamedValue<void> {
public:

  NamedValue( const char * arg_name , T * arg_ptr , size_t arg_size )
    : NamedValue<void>(arg_name,m_ref),
      value(arg_ptr), size( arg_size ),
      m_ref( arg_ptr , arg_size ) {}

  ~NamedValue() {}

  T * const    value ;
  const size_t size ;

private:

  Reference<T*> m_ref ;

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue<T&> : public NamedValue<void> {
public:

  NamedValue( const char * arg_name , T & arg_ref )
    : NamedValue<void>(arg_name,m_ref),
      value(arg_ref), m_ref(arg_ref) {}

  ~NamedValue() {}

  T & value ;

private:

  Reference<T> m_ref ;

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------

template<typename T>
class NamedValue : public NamedValue<void> {
public:

  explicit NamedValue( const char * arg_name )
    : NamedValue<void>( arg_name , m_ref ),
      value(), m_ref(value) {}

  NamedValue( const char * arg_name , const T & arg_val )
    : NamedValue<void>( arg_name , m_ref ),
      value(arg_val), m_ref(value) {}

  ~NamedValue() {}

  T value ;

private:

  Reference<T> m_ref ;

  NamedValue();
  NamedValue( const NamedValue & );
  NamedValue & operator = ( const NamedValue & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace phdmesh

#endif


