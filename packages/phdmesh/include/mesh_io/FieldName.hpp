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

#ifndef phdmesh_FieldName_hpp
#define phdmesh_FieldName_hpp

//----------------------------------------------------------------------

#include <string>

#include <util/CSet.hpp>
#include <mesh/Types.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

class FieldName : public CSetMember<FieldName> {
public:

  /** Encode a descriptive text label given a name, dimension, and indices */
  virtual std::string encode( const std::string & ,
                              const unsigned ,
                              const unsigned * const ,
                              const unsigned * const ) const ;

  /** Decode name and indices from a descriptive text name.
   *  Return if the text matches the descriptor.
   */
  virtual bool decode( const std::string & ,
                             std::string & ,
                       const unsigned ,
                             unsigned * const ) const ;

  virtual const char * name() const ;

  virtual ~FieldName();

  FieldName();

protected:

  void verify( const unsigned ,
               const unsigned * const ,
               const unsigned * const ) const ;

private:
  FieldName( const FieldName & );
  FieldName & operator = ( const FieldName & );
};

const FieldName & array_name();

const FieldName & cartesian_vector();

const FieldName & cylindrical_vector();

void declare( Field<void,0> & , const FieldName & );

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

