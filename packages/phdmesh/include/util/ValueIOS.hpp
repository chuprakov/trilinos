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
 *  class ValueIOS<T> : public ValueIOS<void> {
 *  public:
 *    typedef T ValueType ;
 *    ~ValueIOS();
 *    ValueIOS();
 *    static const ValueIOS<ValueType> & singleton();
 *    virtual void tell(std::ostream &, unsigned, const ValueType &) const ;
 *    virtual void put( std::ostream &, unsigned, const ValueType &) const ;
 *    virtual void get( std::istream &, ValueType & ) const ;
 *  private:
 *    void tellp(std::ostream &, unsigned, const void *) const ;
 *    void putp( std::ostream &, unsigned, const void *) const ;
 *    void getp( std::istream &, void *) const ;
 *  };
 *
 *  The second 'unsigned' argument to the 'tell' and 'put' methods
 *  is a requested indention (number of spaces) to follow any
 *  end-of-line output to the stream.
 *
 *  Specializations are defined here for the following types.
 *    1)  instrinsic_type
 *    2)  instrinsic_type[N]
 *    3)  std::vector< instrinsic_type >
 */
template<class Type> class ValueIOS ;

//----------------------------------------------------------------------

template<>
class ValueIOS<void> {
public:
  virtual ~ValueIOS() {}
  virtual void getp(  std::istream &,                 void *) const = 0 ;
  virtual void putp(  std::ostream &, unsigned, const void *) const = 0 ;
  virtual void tellp( std::ostream &, unsigned, const void *) const = 0 ;

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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Stream io partial specializations for simple scalar arrays

bool value_ios_get_token( std::istream &, int, bool required = false );
bool value_ios_get_block_begin( std::istream & , bool required = false );
bool value_ios_get_block_end(   std::istream & , bool required = false );

void value_ios_get_array(  std::istream & , unsigned , short * );
void value_ios_put_array(  std::ostream & , unsigned , const short * );
void value_ios_tell_array( std::ostream & , unsigned , const short * );

void value_ios_get_array(  std::istream & , unsigned , unsigned short * );
void value_ios_put_array(  std::ostream & , unsigned , const unsigned short * );
void value_ios_tell_array( std::ostream & , unsigned , const unsigned short * );

void value_ios_get_array(  std::istream & , unsigned , int * );
void value_ios_put_array(  std::ostream & , unsigned , const int * );
void value_ios_tell_array( std::ostream & , unsigned , const int * );

void value_ios_get_array(  std::istream & , unsigned , unsigned * );
void value_ios_put_array(  std::ostream & , unsigned , const unsigned * );
void value_ios_tell_array( std::ostream & , unsigned , const unsigned * );

void value_ios_get_array(  std::istream & , unsigned , long * );
void value_ios_put_array(  std::ostream & , unsigned , const long * );
void value_ios_tell_array( std::ostream & , unsigned , const long * );

void value_ios_get_array(  std::istream & , unsigned , unsigned long * );
void value_ios_put_array(  std::ostream & , unsigned , const unsigned long * );
void value_ios_tell_array( std::ostream & , unsigned , const unsigned long * );

void value_ios_get_array(  std::istream & , unsigned , float * );
void value_ios_put_array(  std::ostream & , unsigned , const float * );
void value_ios_tell_array( std::ostream & , unsigned , const float * );

void value_ios_get_array(  std::istream & , unsigned , double * );
void value_ios_put_array(  std::ostream & , unsigned , const double * );
void value_ios_tell_array( std::ostream & , unsigned , const double * );


template<typename T, unsigned N>
class ValueIOS<T[N]> : public ValueIOS<void> {
public:
  typedef T          ScalarType ;
  typedef ScalarType ValueType[N] ;

  virtual void get(  std::istream & s, ValueType & v ) const
    { value_ios_get_array( s , N , & v[0] ); }

  virtual void put(  std::ostream & s, unsigned, const ValueType & v) const
    { value_ios_put_array( s , N , & v[0] ); }

  virtual void tell( std::ostream & s, unsigned, const ValueType & v) const
    { value_ios_tell_array( s , N , & v[0] ); }

  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();

private:

  virtual void tellp( std::ostream & s , unsigned i, const void * p ) const
    { tell( s, i, * reinterpret_cast<const ValueType*>(p) ); }

  virtual void putp(  std::ostream & s , unsigned i, const void * p ) const
    { put( s, i, * reinterpret_cast<const ValueType*>(p) ); }

  virtual void getp(  std::istream & s , void * p ) const
    { get( s, * reinterpret_cast<ValueType*>(p) ); }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Stream io specializations for std::vectors of simple scalar types.

template<>
struct ValueIOS< std::vector<short> > : public ValueIOS<void> {
public:
  typedef std::vector<short> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<unsigned short> > : public ValueIOS<void> {
public:
  typedef std::vector<unsigned short> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<int> > : public ValueIOS<void> {
public:
  typedef std::vector<int> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<unsigned> > : public ValueIOS<void> {
public:
  typedef std::vector<unsigned> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<long> > : public ValueIOS<void> {
public:
  typedef std::vector<long> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<unsigned long> > : public ValueIOS<void> {
public:
  typedef std::vector<unsigned long> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<float> > : public ValueIOS<void> {
public:
  typedef std::vector<float> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

template<>
struct ValueIOS< std::vector<double> > : public ValueIOS<void> {
public:
  typedef std::vector<double> ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
private:
  virtual void tellp( std::ostream & , unsigned , const void * ) const ;
  virtual void putp(  std::ostream & , unsigned , const void * ) const ;
  virtual void getp(  std::istream & , void * ) const ;
};

//----------------------------------------------------------------------

template<>
struct ValueIOS< std::string > : public ValueIOS<void> {
public:
  typedef std::string ValueType ;
  virtual void get(  std::istream &, ValueType & ) const ;
  virtual void put(  std::ostream &, unsigned, const ValueType &) const ;
  virtual void tell( std::ostream &, unsigned, const ValueType &) const ;
  ~ValueIOS() {}
  ValueIOS() : ValueIOS<void>() {}
  static const ValueIOS<ValueType> & singleton();
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
  ~ValueIOS_Quoted() {}
  ValueIOS_Quoted() : ValueIOS<ValueType>() {}
  static const ValueIOS<ValueType> & singleton();
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template<typename T, unsigned N>
const ValueIOS<T[N]> & value_ios_array_singleton()
{ static const ValueIOS<T[N]> io ; return io ; }

}

template<typename T, unsigned N>
inline
const ValueIOS<T[N]> & ValueIOS<T[N]>::singleton()
{ return value_ios_array_singleton<T,N>(); }

}

//----------------------------------------------------------------------

#endif


