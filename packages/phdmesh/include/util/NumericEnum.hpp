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

#ifndef util_NumericEnum_hpp
#define util_NumericEnum_hpp

#include <complex>
#include <util/TypeList.hpp>

namespace phdmesh {

/** List of numeric types, with 'void' as undefined */

typedef TypeList<          void ,
        TypeList< signed   char ,
        TypeList< unsigned char ,
        TypeList< signed   short ,
        TypeList< unsigned short ,
        TypeList< signed   int ,
        TypeList< unsigned int ,
        TypeList< signed   long ,
        TypeList< unsigned long ,
        TypeList<          float ,
        TypeList<          double ,
        TypeList<          std::complex<float> ,
        TypeList<          std::complex<double> ,
        TypeListEnd > > > > > > > > > > > > > NumericTypeList ;

template<typename Type = void> struct NumericEnum ;

template<>
struct NumericEnum<void> {
  enum { OK = StaticAssert< TypeListUnique<NumericTypeList>::value >::OK };

  enum { length  = TypeListLength<NumericTypeList>::value };
  enum { minimum = 1 };
  enum { maximum = length - 1 };

  static const char * name( unsigned ordinal );
  static unsigned     size( unsigned ordinal );

  enum { value = 0 };
};

template<typename Type>
struct NumericEnum {
  enum { value = TypeListIndex< NumericTypeList , Type>::value };

  enum { OK = StaticAssert<
               0 < (int) value &&
                   (int) value < (int) NumericEnum<void>::length >::OK };
};

template<unsigned Ordinal>
struct NumericType {
private:
  enum { OK = StaticAssert< Ordinal < NumericEnum<>::length >::OK };
public:
  typedef typename TypeListAt< NumericTypeList , Ordinal >::type type ;
};


}

#endif

