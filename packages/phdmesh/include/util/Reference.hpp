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

#ifndef util_Reference_hpp
#define util_Reference_hpp

#include <typeinfo>
#include <iosfwd>
#include <limits>
#include <vector>

#include <util/ValueIO.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** @class Reference
 *  @brief Castable reference to a value, array, or vector.
 */
template< typename T = void ,
          template <typename> class IO = ValueIO > class Reference ;

//----------------------------------------------------------------------

template<>
class Reference<void,ValueIO> {
public:

  const std::type_info & type ;
  const size_t           binary_size ;

  //----------------------------------

  virtual       size_t get_max() const = 0 ;
  virtual       size_t put_max() const = 0 ;
  virtual       void * put_void( size_t ) const = 0 ;
  virtual const void * get_void( size_t ) const = 0 ;

  virtual void   write( std::ostream & ) const = 0 ;
  virtual size_t read(  std::istream & ) const = 0 ;
  virtual void   tell(  std::ostream & ) const = 0 ;

  virtual ~Reference();

  //----------------------------------

  template<typename T>
  void assert_type() const
    { if ( type != typeid(T) ) { throw_type( typeid(T) ); } }

  template<typename T>
  const T & get( size_t i = 0 ) const
    { assert_type<T>(); return * reinterpret_cast<const T*>( get_void(i) ); }

  template<typename T>
  T & put( size_t i = 0 ) const
    { assert_type<T>(); return * reinterpret_cast<T*>( put_void(i) ); }

protected:

  Reference( const std::type_info & t , size_t s ) : type(t), binary_size(s) {}

private:

  void throw_type( const std::type_info & ) const ;

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};

//----------------------------------------------------------------------

template<typename T>
class Reference< const std::vector<T> , ValueIO > : public Reference<> {
public:
  const std::vector<T> & value ;

  explicit Reference( const std::vector<T> & arg )
    : Reference<>( typeid(T) , ValueIO<T>::binary_size ), value(arg) {}

  ~Reference() {}

  void write( std::ostream & s ) const
    { ValueIO<T>::write( s , & value[0] , value.size() ); }

  size_t read( std::istream & ) const { return 0 ; }

  void tell( std::ostream & s ) const
    { ValueIO<T>::tell( s , & value[0] , get_max() , put_max() ); }

private:

  size_t get_max() const { return value.size(); }
  size_t put_max() const { return 0 ; }

  void * put_void( size_t ) const { return NULL ; }

  const void * get_void( size_t i ) const
    { return i < value.size() ? & value[i] : NULL ; }

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};


template<typename T>
class Reference< std::vector<T> , ValueIO > : public Reference<> {
public:
  std::vector<T> & value ;

  explicit Reference( std::vector<T> & arg )
    : Reference<>( typeid(T) , ValueIO<T>::binary_size ), value(arg) {}

  ~Reference() {}

  void write( std::ostream & s ) const
    { ValueIO<T>::write( s , & value[0] , value.size() ); }

  size_t read( std::istream & s ) const
    {
      size_t i = ValueIO<T>::read( s , & value[0] , value.size() );
      if ( i == value.size() ) {
        T tmp ;
        while ( ValueIO<T>::read(s, &tmp, 1) )
          { value.push_back(tmp); ++i ; }
      }
      return i ;
    }

  void tell( std::ostream & s ) const
    { ValueIO<T>::tell( s , & value[0] , get_max() , put_max() ); }

private:

  size_t get_max() const { return value.size(); }
  size_t put_max() const { return std::numeric_limits<size_t>::max(); }

  void * put_void( size_t i ) const
    {
      if ( value.size() <= i ) { value.resize(i+1); }
      return i < value.size() ? & value[i] : NULL ;
    }

  const void * get_void( size_t i ) const
    { return i < value.size() ? & value[i] : NULL ; }

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};

//----------------------------------------------------------------------

template<typename T>
class Reference<const T*> : public Reference<> {
public:
  const T * const value ;
  const size_t    size ;

  Reference( const T * arg_ptr , size_t arg_size )
    : Reference<>( typeid(T) , ValueIO<T>::binary_size ),
      value(arg_ptr), size(arg_size) {}

  ~Reference() {}

  void write( std::ostream & s ) const
    { ValueIO<T>::write( s , value , size ); }

  size_t read( std::istream & ) const { return 0 ; }

  void tell( std::ostream & s ) const
    { ValueIO<T>::tell( s , value , get_max() , put_max() ); }

private:

        size_t get_max() const { return size ; }
        size_t put_max() const { return 0 ; }
        void * put_void(size_t )  const {return NULL ; }
  const void * get_void(size_t i) const {return i < size ? value+i : NULL;}

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};


template<typename T>
class Reference<T*> : public Reference<> {
public:
  T * const    value ;
  const size_t size ;

  Reference( T * arg_ptr , size_t arg_size )
    : Reference<>( typeid(T) , ValueIO<T>::binary_size ),
      value(arg_ptr), size( arg_size ) {}

  ~Reference() {}

  void write( std::ostream & s ) const
    { ValueIO<T>::write( s , value , size ); }

  size_t read( std::istream & s ) const
    { return ValueIO<T>::read( s , value , size ); }

  void tell( std::ostream & s ) const
    { ValueIO<T>::tell( s , value , get_max() , put_max() ); }

private:

        size_t get_max() const { return size ; }
        size_t put_max() const { return size ; }
        void * put_void(size_t i) const {return i < size ? value+i : NULL;}
  const void * get_void(size_t i) const {return i < size ? value+i : NULL;}

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};

//----------------------------------------------------------------------

template<typename T>
class Reference<const T> : public Reference<> {
public:
  const T & value ;

  explicit Reference( const T & arg )
    : Reference<>( typeid(T) , ValueIO<T>::binary_size ), value(arg) {}

  ~Reference() {}

  void write( std::ostream & s ) const
    { ValueIO<T>::write( s , & value , 1 ); }

  size_t read( std::istream & ) const { return 0 ; }

  void tell( std::ostream & s ) const
    { ValueIO<T>::tell( s , & value , get_max() , put_max() ); }

private:

        size_t get_max() const { return 1 ; }
        size_t put_max() const { return 0 ; }
        void * put_void( size_t i ) const { return NULL ; }
  const void * get_void( size_t i ) const { return i ? NULL : & value ; }

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};

template<typename T, template <typename> class IO >
class Reference : public Reference<> {
public:
  T & value ;

  explicit Reference( T & arg )
    : Reference<>( typeid(T) , IO<T>::binary_size ), value(arg) {}

  ~Reference() {}

  void write( std::ostream & s ) const
    { IO<T>::write( s , & value , 1 ); }

  size_t read( std::istream & s ) const
    { return IO<T>::read( s , & value , 1 ); }

  void tell( std::ostream & s ) const
    { IO<T>::tell( s , & value , get_max() , put_max() ); }

private:

        size_t get_max() const { return 1 ; }
        size_t put_max() const { return 1 ; }
        void * put_void( size_t i ) const { return i ? NULL : & value ; }
  const void * get_void( size_t i ) const { return i ? NULL : & value ; }

  Reference();
  Reference( const Reference & );
  Reference & operator = ( const Reference & );
};

//----------------------------------------------------------------------

}

#endif

