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

#ifndef util_ValueIO_hpp
#define util_ValueIO_hpp

#include <iosfwd>

namespace phdmesh {

//----------------------------------------------------------------------
/** @class  ValueIO
 *  @brief  Type specific bundle of stream input, output, and
 *          description functions.
 *
 *  The ValueIOS template class is intended to bundle stream
 *  input, output, and description methods for the templated type.
 */
template<typename T> class ValueIO ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Stream io specializations for simple scalar values

template<>
class ValueIO<unsigned char> {
public:
  typedef unsigned char ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<signed char> {
public:
  typedef signed char ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<short> {
public:
  typedef short ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<unsigned short> {
public:
  typedef unsigned short ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<int> {
public:
  typedef int ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<unsigned int> {
public:
  typedef unsigned int ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<long> {
public:
  typedef long ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<unsigned long> {
public:
  typedef unsigned long ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<float> {
public:
  typedef float ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

template<>
class ValueIO<double> {
public:
  typedef double ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Stream io specialization support for enumerated values
// Functions accept array of ValueIO_Enum with final member
// having a 'name' value of NULL.  Name comparison is
// case-insensitive.

std::string read_name( std::istream & );

struct ValueIO_Enum {
  const char * name ;
        long   value ;
};

const ValueIO_Enum * find( const ValueIO_Enum * , const char * );
const ValueIO_Enum * find( const ValueIO_Enum * , long );
const ValueIO_Enum * read( const ValueIO_Enum * , std::istream & );
void                 tell( const ValueIO_Enum * , std::ostream & );

template<>
class ValueIO<bool> {
public:
  typedef bool ValueType ;
  enum { binary_size = sizeof(ValueType) };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/* The usual non-white space string */

template<>
struct ValueIO< std::string > {
public:
  typedef std::string ValueType ;
  enum { binary_size = 0 /* Cannot perform simple binary IO */ };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

/* A string containing a name pattern */

template <typename T> struct ValueIO_Name ;

template<>
struct ValueIO_Name<std::string> {
public:
  typedef std::string ValueType ;
  enum { binary_size = 0 /* Cannot perform simple binary IO */ };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

/* A string enclosed in quotes */

template <typename T> struct ValueIO_Quoted ;

template<>
struct ValueIO_Quoted<std::string> {
public:
  typedef std::string ValueType ;
  enum { binary_size = 0 /* Cannot perform simple binary IO */ };
  static size_t read(  std::istream &,       ValueType * , size_t );
  static void  write( std::ostream &, const ValueType * , size_t );
  static void  tell(  std::ostream &, const ValueType * , size_t , size_t );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

#endif


