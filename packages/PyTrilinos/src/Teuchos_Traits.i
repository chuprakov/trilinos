// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"

// Specializations of Teuchos::SerializationTraits<>:
//     unsigned char, short, unsigned short, unsigned int, long,
//     unsigned long, long long, unsigned long long
namespace Teuchos
{
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned char>
  : public DirectSerializationTraits<Ordinal,unsigned char>
{};
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned short>
  : public DirectSerializationTraits<Ordinal,unsigned short>
{};
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned int>
  : public DirectSerializationTraits<Ordinal,unsigned int>
{};
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned long>
  : public DirectSerializationTraits<Ordinal,unsigned long>
{};
template<typename Ordinal>
class SerializationTraits<Ordinal,long long>
  : public DirectSerializationTraits<Ordinal,long long>
{};
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned long long>
  : public DirectSerializationTraits<Ordinal,unsigned long long>
{};
}
%}

///////////////////////////////////
// Teuchos::ScalarTraits support //
///////////////////////////////////
%include "Teuchos_ScalarTraits.hpp"
// Specialization of ScalarTraits<> for type long
%inline 
%{
namespace Teuchos
{
// Type unsigned char
template<>
struct ScalarTraits<unsigned char>
{
  typedef unsigned char magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned char a)
  { return static_cast<unsigned char>(std::fabs(static_cast<double>(a))); }
  static inline unsigned char zero()  { return 0; }
  static inline unsigned char one()   { return 1; }
  static inline unsigned char conjugate(unsigned char x) { return x; }
  static inline unsigned char real(unsigned char x) { return x; }
  static inline unsigned char imag(unsigned char x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned char random() { return std::rand(); }
  static inline std::string name() { return "unsigned char"; }
  static inline unsigned char squareroot(unsigned char x)
  { return (unsigned char) std::sqrt((double) x); }
  static inline unsigned char pow(unsigned char x, unsigned char y)
  { return (unsigned char) std::pow((double)x,(double)y); }
};

// Type unsigned short
template<>
struct ScalarTraits<unsigned short>
{
  typedef unsigned short magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned short a)
  { return static_cast<unsigned short>(std::fabs(static_cast<double>(a))); }
  static inline unsigned short zero()  { return 0; }
  static inline unsigned short one()   { return 1; }
  static inline unsigned short conjugate(unsigned short x) { return x; }
  static inline unsigned short real(unsigned short x) { return x; }
  static inline unsigned short imag(unsigned short x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned short random() { return std::rand(); }
  static inline std::string name() { return "unsigned short"; }
  static inline unsigned short squareroot(unsigned short x)
  { return (unsigned short) std::sqrt((double) x); }
  static inline unsigned short pow(unsigned short x, unsigned short y)
  { return (unsigned short) std::pow((double)x,(double)y); }
};

// Type unsigned int
template<>
struct ScalarTraits<unsigned int>
{
  typedef unsigned int magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned int a)
  { return static_cast<unsigned int>(std::fabs(static_cast<double>(a))); }
  static inline unsigned int zero()  { return 0; }
  static inline unsigned int one()   { return 1; }
  static inline unsigned int conjugate(unsigned int x) { return x; }
  static inline unsigned int real(unsigned int x) { return x; }
  static inline unsigned int imag(unsigned int x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned int random() { return std::rand(); }
  static inline std::string name() { return "unsigned int"; }
  static inline unsigned int squareroot(unsigned int x)
  { return (unsigned int) std::sqrt((double) x); }
  static inline unsigned int pow(unsigned int x, unsigned int y)
  { return (unsigned int) std::pow((double)x,(double)y); }
};

// Type unsigned long
template<>
struct ScalarTraits<unsigned long>
{
  typedef unsigned long magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned long a)
  { return static_cast<unsigned long>(std::fabs(static_cast<double>(a))); }
  static inline unsigned long zero()  { return 0; }
  static inline unsigned long one()   { return 1; }
  static inline unsigned long conjugate(unsigned long x) { return x; }
  static inline unsigned long real(unsigned long x) { return x; }
  static inline unsigned long imag(unsigned long x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned long random() { return std::rand(); }
  static inline std::string name() { return "unsigned long"; }
  static inline unsigned long squareroot(unsigned long x)
  { return (unsigned long) std::sqrt((double) x); }
  static inline unsigned long pow(unsigned long x, unsigned long y)
  { return (unsigned long) std::pow((double)x,(double)y); }
};

// Type long long
template<>
struct ScalarTraits<long long>
{
  typedef long long magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(long long a)
  { return static_cast<long long>(std::fabs(static_cast<double>(a))); }
  static inline long long zero()  { return 0; }
  static inline long long one()   { return 1; }
  static inline long long conjugate(long long x) { return x; }
  static inline long long real(long long x) { return x; }
  static inline long long imag(long long x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline long long random() { return (-1 + 2*rand()); }
  static inline std::string name() { return "long long"; }
  static inline long long squareroot(long long x)
  { return (long long) std::sqrt((double) x); }
  static inline long long pow(long long x, long long y)
  { return (long long) std::pow((double)x,(double)y); }
};

// Type unsigned long long
template<>
struct ScalarTraits<unsigned long long>
{
  typedef unsigned long long magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned long long a)
  { return static_cast<unsigned long long>(std::fabs(static_cast<double>(a))); }
  static inline unsigned long long zero()  { return 0; }
  static inline unsigned long long one()   { return 1; }
  static inline unsigned long long conjugate(unsigned long long x) { return x; }
  static inline unsigned long long real(unsigned long long x) { return x; }
  static inline unsigned long long imag(unsigned long long x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned long long random() { return std::rand(); }
  static inline std::string name() { return "unsigned long long"; }
  static inline unsigned long long squareroot(unsigned long long x)
  { return (unsigned long long) std::sqrt((double) x); }
  static inline unsigned long long pow(unsigned long long x, unsigned long long y)
  { return (unsigned long long) std::pow((double)x,(double)y); }
};

}
%}
%template(ScalarTraitsByte     ) Teuchos::ScalarTraits< char               >;
%template(ScalarTraitsUbyte    ) Teuchos::ScalarTraits< unsigned char      >;
%template(ScalarTraitsShort    ) Teuchos::ScalarTraits< short              >;
%template(ScalarTraitsUshort   ) Teuchos::ScalarTraits< unsigned short     >;
%template(ScalarTraitsInt      ) Teuchos::ScalarTraits< int                >;
%template(ScalarTraitsUint     ) Teuchos::ScalarTraits< unsigned int       >;
%template(ScalarTraitsLong     ) Teuchos::ScalarTraits< long               >;
%template(ScalarTraitsUlong    ) Teuchos::ScalarTraits< unsigned long      >;
%template(ScalarTraitsLongLong ) Teuchos::ScalarTraits< long long          >;
%template(ScalarTraitsUlongLong) Teuchos::ScalarTraits< unsigned long long >;
%template(ScalarTraitsFloat    ) Teuchos::ScalarTraits< float              >;
%template(ScalarTraitsDouble   ) Teuchos::ScalarTraits< double             >;
%pythoncode
%{
def ScalarTraits(scalarType):
    """
    ScalarTraits(str scalarType) -> ScalarTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which traits are requested.  Limited NumPy-style type
    specification is supported:

        'b'  byte
        'B'  unsigned byte
        'h'  short
        'H'  unsigned short
        'i'  int
        'I'  unsigned int
        'l'  long
        'L'  unsigned long
        'q'  long long
        'Q'  unsigned long long
        'f'  float
        'd'  double
    """
    if   scalarType == 'b': return ScalarTraitsByte()
    elif scalarType == 'B': return ScalarTraitsUbyte()
    elif scalarType == 'h': return ScalarTraitsShort()
    elif scalarType == 'H': return ScalarTraitsUshort()
    elif scalarType == 'i': return ScalarTraitsInt()
    elif scalarType == 'I': return ScalarTraitsUint()
    elif scalarType == 'l': return ScalarTraitsLong()
    elif scalarType == 'L': return ScalarTraitsUlong()
    elif scalarType == 'q': return ScalarTraitsLongLong()
    elif scalarType == 'Q': return ScalarTraitsUlongLong()
    elif scalarType == 'f': return ScalarTraitsFloat()
    elif scalarType == 'd': return ScalarTraitsDouble()
    raise NotImplementedError, "ScalarTraits for " + repr(scalarType) + " not supported"
%}

//////////////////////////////////////////
// Teuchos::SerializationTraits support //
//////////////////////////////////////////
%define
%teuchos_serialtraits_for_type(type,Type)
%ignore Teuchos::DirectSerializationTraits<long,type>::convertToCharPtr(const type*);
%ignore Teuchos::DirectSerializationTraits<long,type>::convertToCharPtr(const char*);
%ignore Teuchos::DirectSerializationTraits<long,type>::convertFromCharPtr(const char*);
%template(DirectSerializationTraits ## Type)
    Teuchos::DirectSerializationTraits<long, type>;
%ignore Teuchos::SerializationTraits<long,type>::convertToCharPtr(const type*);
%ignore Teuchos::SerializationTraits<long,type>::convertToCharPtr(const char*);
%ignore Teuchos::SerializationTraits<long,type>::convertFromCharPtr(const char*);
%template(SerializationTraits ## Type) Teuchos::SerializationTraits<long, type>;
%enddef
%include "Teuchos_SerializationTraits.hpp"
%teuchos_serialtraits_for_type(char              , Byte     )
%teuchos_serialtraits_for_type(unsigned char     , Ubyte    )
%teuchos_serialtraits_for_type(short             , Short    )
%teuchos_serialtraits_for_type(unsigned short    , Ushort   )
%teuchos_serialtraits_for_type(int               , Int      )
%teuchos_serialtraits_for_type(unsigned int      , Uint     )
%teuchos_serialtraits_for_type(long              , Long     )
%teuchos_serialtraits_for_type(unsigned long     , Ulong    )
%teuchos_serialtraits_for_type(long long         , LongLong )
%teuchos_serialtraits_for_type(unsigned long long, UlongLong)
%teuchos_serialtraits_for_type(float             , Float    )
%teuchos_serialtraits_for_type(double            , Double   )
%pythoncode
%{
def DirectSerializationTraits(scalarType):
    """
    DirectSerializationTraits(str scalarType) -> DirectSerializationTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which direct serialization is requested.  Limited NumPy-style type
    specification is supported:

        'b'  byte
        'B'  unsigned byte
        'h'  short
        'H'  unsigned short
        'i'  int
        'I'  unsigned int
        'l'  long
        'L'  unsigned long
        'q'  long long
        'Q'  unsigned long long
        'f'  float
        'd'  double
    """
    if   scalarType == 'b': return DirectSerializationTraitsByte()
    elif scalarType == 'B': return DirectSerializationTraitsUbyte()
    elif scalarType == 'h': return DirectSerializationTraitsShort()
    elif scalarType == 'H': return DirectSerializationTraitsUshort()
    elif scalarType == 'i': return DirectSerializationTraitsInt()
    elif scalarType == 'I': return DirectSerializationTraitsUint()
    elif scalarType == 'l': return DirectSerializationTraitsLong()
    elif scalarType == 'L': return DirectSerializationTraitsUlong()
    elif scalarType == 'q': return DirectSerializationTraitsLongLong()
    elif scalarType == 'Q': return DirectSerializationTraitsUlongLong()
    elif scalarType == 'f': return DirectSerializationTraitsFloat()
    elif scalarType == 'd': return DirectSerializationTraitsDouble()
    raise NotImplementedError, "DirectSerializationTraits for " + repr(scalarType) + \
	  " not supported"

def SerializationTraits(scalarType):
    """
    SerializationTraits(str scalarType) -> SerializationTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which serialization is requested.  Limited NumPy-style type
    specification is supported:

        'b'  byte
        'B'  unsigned byte
        'h'  short
        'H'  unsigned short
        'i'  int
        'I'  unsigned int
        'l'  long
        'L'  unsigned long
        'q'  long long
        'Q'  unsigned long long
        'f'  float
        'd'  double
    """
    if   scalarType == 'b': return SerializationTraitsByte()
    elif scalarType == 'B': return SerializationTraitsUbyte()
    elif scalarType == 'h': return SerializationTraitsShort()
    elif scalarType == 'H': return SerializationTraitsUshort()
    elif scalarType == 'i': return SerializationTraitsInt()
    elif scalarType == 'I': return SerializationTraitsUint()
    elif scalarType == 'l': return SerializationTraitsLong()
    elif scalarType == 'L': return SerializationTraitsUlong()
    elif scalarType == 'q': return SerializationTraitsLongLong()
    elif scalarType == 'Q': return SerializationTraitsUlongLong()
    elif scalarType == 'f': return SerializationTraitsFloat()
    elif scalarType == 'd': return SerializationTraitsDouble()
    raise NotImplementedError, "SerializationTraits for " + repr(scalarType) + \
	  " not supported"
%}
