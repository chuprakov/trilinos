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

#ifndef util_ValueIOS_hpp
#define util_ValueIOS_hpp

#include <iosfwd>
#include <typeinfo>
#include <vector>

namespace phdmesh {

//----------------------------------------------------------------------
/** @class  ValueIOS
 *  @brief  Type specific bundle of stream input, output, and
 *          description functions.
 *
 *  The ValueIOS template class is intended to bundle stream
 *  input, output, and description methods for the templated type.
 *  Two version of these methods are required: (1) type specific
 *  and (2) type anonymous.  The type anonymous versions enable
 *  heterogeneously typed containers (e.g. NamedValueSet) to
 *  perform stream io on their members.
 *  
 *  The required template class specialization for a type as follows.
 *
 *  template<>
 *  class ValueIOS<ScalarType> : public ValueIOS<void> {
 *  public:
 *    typedef ScalarType ValueType ;
 *    ~ValueIOS();
 *    ValueIOS();
 *    virtual void tell(std::ostream &, unsigned, const ValueType &) const ;
 *    virtual void put( std::ostream &, unsigned, const ValueType &) const ;
 *    virtual void get( std::istream &, ValueType & ) const ;
 *    virtual const std::type_info & type() const ;
 *  private:
 *    void tellp(std::ostream &, unsigned, const void *) const ;
 *    void putp( std::ostream &, unsigned, const void *) const ;
 *    void getp( std::istream &, void *) const ;
 *  };
 *
 *  The second 'unsigned' argument to the 'tell' and 'put' methods
 *  is a requested indention (number of spaces) to follow any
 *  end-of-line output to the stream.
 */
template<class ScalarType> class ValueIOS ;

class ValueIOSPolicy ;

//----------------------------------------------------------------------
/** @class ValueIOSPolicy
 *  @brief A collection of stream io operations for scalar types.
 */
class ValueIOSPolicy {
public:

  /** Set scalar value stream io for the given scalar type */
  void replace( const ValueIOS<void> & );

  /** Get scalar value stream io for the given scalar type */
  const ValueIOS<void> * get_void( const std::type_info & ) const ;

  /** Get scalar value stream io for the given scalar type */
  template<class ScalarType>
  const ValueIOS<ScalarType> * get() const
    {
      return static_cast<const ValueIOS<ScalarType>*>(
        get_void( typeid(ScalarType) ) );
    }

  ~ValueIOSPolicy();
  ValueIOSPolicy();
  ValueIOSPolicy( const ValueIOSPolicy & );

private:

  ValueIOSPolicy & operator = ( const ValueIOSPolicy & );

  std::vector<const ValueIOS<void> *> m_ios ;
};

//----------------------------------------------------------------------

template<>
class ValueIOS<void> {
public:
  virtual ~ValueIOS() {}
  virtual void getp(  std::istream &,                 void *) const = 0 ;
  virtual void putp(  std::ostream &, unsigned, const void *) const = 0 ;
  virtual void tellp( std::ostream &, unsigned, const void *) const = 0 ;
  virtual const std::type_info & type() const = 0 ;
protected:
  ValueIOS() {}
private:
  ValueIOS( const ValueIOS<void> & );
  ValueIOS<void> & operator = ( const ValueIOS<void> & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Stream io specializations for simple scalar values

template<>
class ValueIOS<short> : public ValueIOS<void> {
public:
  typedef short ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<unsigned short> : public ValueIOS<void> {
public:
  typedef unsigned short ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<int> : public ValueIOS<void> {
public:
  typedef int ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<unsigned int> : public ValueIOS<void> {
public:
  typedef unsigned int ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<long> : public ValueIOS<void> {
public:
  typedef long ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<unsigned long> : public ValueIOS<void> {
public:
  typedef unsigned long ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<float> : public ValueIOS<void> {
public:
  typedef float ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
class ValueIOS<double> : public ValueIOS<void> {
public:
  typedef double ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Stream io specialization support for enumerated values
// Functions accept array of ValueIOS_Enum with final member
// having a 'm_name' value of NULL.  Name comparison is
// case-insensitive.

struct ValueIOS_Enum {
  const char * m_name ;
        long   m_value ;
};

long         enum_value_of_name( const ValueIOS_Enum * , const char * );
const char * enum_name_of_value( const ValueIOS_Enum * , const long );
void         enum_tell( std::ostream & , unsigned ,
                        const char * , const ValueIOS_Enum * );

template<>
class ValueIOS<bool> : public ValueIOS<void> {
public:
  typedef bool ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
struct ValueIOS< std::string > : public ValueIOS<void> {
public:
  typedef std::string ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

struct ValueIOS_Quoted : public ValueIOS< std::string > {
public:
  typedef ValueIOS< std::string >::ValueType ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  virtual const std::type_info & type() const { return typeid(ValueType); }
  ~ValueIOS_Quoted() {}
  ValueIOS_Quoted() : ValueIOS<ValueType>() {}
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

#endif


