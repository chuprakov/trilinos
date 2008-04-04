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

#ifndef util_Basics_hpp
#define util_Basics_hpp

#include <limits>
#include <phdmesh_config.h>

namespace phdmesh {

//----------------------------------------------------------------------
// Compile time assertion
//     enum { ok = StaticAssert< logical_expression >::OK };
//     StaticAssert< logical_expression >::ok();
//  For logical_expression == true  it generates a valid no-op
//  For logical_expression == false it generates a compile error

template<bool> struct StaticAssert ;

template<> struct StaticAssert<true> {
  enum { OK = true };
  static bool ok() { return true ; }
};

template<> struct StaticAssert<false> {};

//----------------------------------------------------------------------
// Compile time comparison of types

template<typename T1, typename T2> struct SameType ;

template<typename T> struct SameType<T,T> { enum { value = true }; };

template <typename T1, typename T2> struct SameType { enum { value = false }; };

//----------------------------------------------------------------------
// Selection of an integer type based upon sign and size

template<
  unsigned N_char  = std::numeric_limits<unsigned char >::digits ,
  unsigned N_short = std::numeric_limits<unsigned short>::digits ,
  unsigned N_int   = std::numeric_limits<unsigned int  >::digits ,
  unsigned N_long  = std::numeric_limits<unsigned long >::digits >
struct IntegerTypes ;

template<>
struct IntegerTypes<8,16,32,64> {
  typedef signed   char  int8_type ;
  typedef unsigned char  uint8_type ;
  typedef signed   short int16_type ;
  typedef unsigned short uint16_type ;
  typedef signed   int   int32_type ;
  typedef unsigned int   uint32_type ;
  typedef signed   long  int64_type ;
  typedef unsigned long  uint64_type ;
};

template<
  unsigned N_char ,
  unsigned N_short ,
  unsigned N_int ,
  unsigned N_long >
struct IntegerTypes {
  typedef signed   char      int8_type ;
  typedef unsigned char      uint8_type ;
  typedef signed   short     int16_type ;
  typedef unsigned short     uint16_type ;
  typedef signed   int       int32_type ;
  typedef unsigned int       uint32_type ;
  typedef signed   long long int64_type ;
  typedef unsigned long long uint64_type ;

  enum { uint8_digits  = std::numeric_limits<uint8_type >::digits };
  enum { uint16_digits = std::numeric_limits<uint16_type>::digits };
  enum { uint32_digits = std::numeric_limits<uint32_type>::digits };
  enum { uint64_digits = std::numeric_limits<uint64_type>::digits };

  enum { OK_8  = StaticAssert<  8 == uint8_digits  >::OK };
  enum { OK_16 = StaticAssert< 16 == uint16_digits >::OK };
  enum { OK_32 = StaticAssert< 32 == uint32_digits >::OK };
  enum { OK_64 = StaticAssert< 64 == uint64_digits >::OK };
};

typedef IntegerTypes<>::int8_type    int8_type ;
typedef IntegerTypes<>::uint8_type   uint8_type ;
typedef IntegerTypes<>::int16_type   int16_type ;
typedef IntegerTypes<>::uint16_type  uint16_type ;
typedef IntegerTypes<>::int32_type   int32_type ;
typedef IntegerTypes<>::uint32_type  uint32_type ;
typedef IntegerTypes<>::int64_type   int64_type ;
typedef IntegerTypes<>::uint64_type  uint64_type ;

} // namespace phdmesh

#endif

