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

#ifndef phdmesh_ExoII_hpp
#define phdmesh_ExoII_hpp

#include <mesh/Types.hpp>

namespace phdmesh {
namespace exodus {

enum ElementType {
  CIRCLE /* Attributes: radius */ ,
  SPHERE /* Attributes: radius */ ,
  TRUSS  /* Attributes: cross-section area */ ,
  BEAM   /* Attributes: cross-section area and moments */ ,
  SHELL  /* Attributes: thickness */ ,
  QUAD , TRIANGLE , PYRAMID , TETRA , WEDGE , HEX ,
  UNDEFINED };

class FilePart ;

//----------------------------------------------------------------------

class FileSchema {
public:

  ~FileSchema();

  FileSchema( Schema                & arg_schema ,
              const Field<double,1> & arg_node_coordinates ,
              const Field<double,1> & arg_elem_attributes ,
              const unsigned          arg_writer_rank = 0 );

  FileSchema( Schema                & arg_schema ,
              const Field<double,1> & arg_node_coordinates ,
              const Field<double,1> & arg_elem_attributes ,
              const std::string     & arg_file_path ,
              ParallelMachine         arg_comm ,
              const unsigned          arg_writer_rank = 0 );

  /** Declare element part, default number of attributes. */
  Part & declare_part( const std::string & arg_name ,
                       int                 arg_id ,
                       ElementType         arg_element_type ,
                       unsigned            arg_number_nodes ,
                       unsigned            arg_num_attributes = 0 );

  /** Declare a node, edge, or face part */
  Part & declare_part( const std::string & arg_name ,
                       int                 arg_id ,
                       EntityType          arg_type );

  /** Assign contiguous global indices [1..#] to nodes and elements.
   *  Elements are ordered by element block and then by identifier.
   */
  void assign_indices( Mesh & ) const ;

  Schema                & m_schema ;
  const unsigned          m_io_rank ;
  const unsigned          m_dimension ;
  const Field<double,1> & m_field_node_coord ;
  const Field<double,1> & m_field_elem_attr ;
  const Field<int,1>    & m_field_node_index ;
  const Field<int,1>    & m_field_edge_index ;
  const Field<int,1>    & m_field_face_index ;
  const Field<int,1>    & m_field_elem_index ;

  const std::vector<const FilePart*> & parts( EntityType t ) const
    { return m_parts[t] ; }

private:
  FileSchema();
  FileSchema( const FileSchema & );
  FileSchema & operator = ( const FileSchema & );

  std::vector<const FilePart*> m_parts[ EntityTypeMaximum ];
};

//----------------------------------------------------------------------

struct FieldIO {
  const Field<void,0> * m_field ;
  unsigned              m_offset ;
  int                   m_var_index ;
  const FilePart      * m_part ;
};


class FileOutput {
public:
  ~FileOutput();

  /** Create an output file for a collection of fields. */
  FileOutput( const FileSchema & ,
              const Mesh & ,
              const std::string & arg_file_path ,
              const std::string & arg_title ,
              const bool          arg_storage_double ,
              const std::vector< const Field<void,0> * > & ,
              const int * const arg_processor = NULL );

  /** Write a snapshot of field values */
  void write( double );

  const FileSchema & m_schema ;
  const Mesh       & m_mesh ;

  int exo_id() const { return m_exo_id ; }
  int exo_step() const { return m_counter ; }

  const std::vector<int> & global_counts() const ;
  const std::vector< FieldIO > & field_node_universal() const
    { return m_field_node_universal ; }

  const std::vector< FieldIO > & field_elem() const
    { return m_field_elem ; }

private:
  FileOutput();
  FileOutput( const FileOutput & );
  FileOutput & operator = ( const FileOutput & );

  int m_exo_id ;
  int m_counter ;
  int m_max_buffer ;

  std::vector< FieldIO > m_field_node_universal ;
  std::vector< FieldIO > m_field_elem ;
  std::vector<int> m_global_counts ;
};

//----------------------------------------------------------------------

class FileInput {
public:
  ~FileInput();

  FileInput( const FileSchema & , Mesh & ,
             const std::string & arg_file_path ,
             const std::vector< const Field<void,0> * > & );

  double read();

  const FileSchema & m_schema ;
        Mesh       & m_mesh ;

  int exo_id() const { return m_exo_id ; }
  int exo_step() const { return m_counter ; }

private:
  FileInput();
  FileInput( const FileInput & );
  FileInput & operator = ( const FileInput & );

  int m_exo_id ;
  int m_counter ;
  int m_max_buffer ;

  std::vector< FieldIO > m_field_node_universal ;
  std::vector< FieldIO > m_field_elem ;
};

} // namespace exodus
} // namespace phdmesh

//----------------------------------------------------------------------

#endif

