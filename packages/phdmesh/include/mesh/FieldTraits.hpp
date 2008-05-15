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
 * @author H. Carter Edwards
 */

#ifndef phdmesh_FieldTraits_hpp
#define phdmesh_FieldTraits_hpp

//----------------------------------------------------------------------

#include <util/Dimension.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

struct Cartesian : public DimensionTraits {

  const char * name() const ;

  std::string encode( unsigned size , unsigned index ) const ;

  unsigned decode( unsigned size , const std::string & ) const ;

  static const DimensionTraits * descriptor();

private:
  Cartesian() {}
  Cartesian( const Cartesian & );
  Cartesian & operator = ( const Cartesian & );
};

struct Cylindrical : public DimensionTraits {

  const char * name() const ;

  std::string encode( unsigned size , unsigned index ) const ;

  unsigned decode( unsigned size , const std::string & ) const ;

  static const DimensionTraits * descriptor();

private:
  Cylindrical() {}
  Cylindrical( const Cylindrical & );
  Cylindrical & operator = ( const Cylindrical & );
};

}

#endif

