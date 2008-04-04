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

#include <string.h>
#include <stdio.h>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <util/ParallelReduce.hpp>
#include <util/ParallelComm.hpp>

#include <mesh/EntityType.hpp>
#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>

#include <mesh_io/ExoII.hpp>
#include <mesh_io/FieldName.hpp>

//----------------------------------------------------------------------

#if defined( PHDMESH_HAS_MPI )

#include <mpi.h>

namespace {

void allgather( phdmesh::ParallelMachine c ,
                unsigned * local , unsigned * global , int n )
{
  MPI_Allgather( local , n , MPI_UNSIGNED , global , n , MPI_UNSIGNED , c );
}

}

#else

namespace {

void allgather( phdmesh::ParallelMachine ,
                unsigned * local , unsigned * global , int n )
{
  for ( int i = 0 ; i < n ; ++i ) { global[i] = local[i] ; }
}

}

#endif

//----------------------------------------------------------------------

namespace phdmesh {
namespace exodus {

//----------------------------------------------------------------------

class FilePart {
public:
  Part            & m_part ;
  const int         m_identifier ;
  const EntityType  m_entity_type ;
  const ElementType m_element_type ;
  const int         m_number_nodes ;
  const int         m_number_attr ;

  FilePart( Part      & arg_part ,
            int         arg_id ,
            EntityType  arg_entity_type ,
            ElementType arg_element_type ,
            int         arg_number_nodes ,
            int         arg_number_attributes )
    : m_part( arg_part ),
      m_identifier( arg_id ),
      m_entity_type( arg_entity_type ),
      m_element_type( arg_element_type ),
      m_number_nodes( arg_number_nodes ),
      m_number_attr( arg_number_attributes )
    {}
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// Exodus indices: [ global , local ]
const Field<int,1> & exo_index( Schema & arg_schema , EntityType arg_type )
{
  static const char name[] = "exodusII_local_index" ;
  Field<int,1> & f =
    arg_schema.declare_field<int,1>( arg_type , std::string(name) , 1 );
  arg_schema.declare_field_dimension( f , arg_schema.universal_part() , 2 );
  return f ;
}

}

FileSchema::FileSchema(
  Schema                & arg_schema ,
  const Field<double,1> & arg_node_coordinates ,
  const Field<double,1> & arg_elem_attributes ,
  const unsigned          arg_writer_rank )
  : m_schema( arg_schema ),
    m_io_rank( arg_writer_rank ),
    m_dimension( arg_node_coordinates.max_length() ),
    m_field_node_coord( arg_node_coordinates ),
    m_field_elem_attr(  arg_elem_attributes ),
    m_field_node_index( exo_index( arg_schema , Node ) ),
    m_field_edge_index( exo_index( arg_schema , Edge ) ),
    m_field_face_index( exo_index( arg_schema , Face ) ),
    m_field_elem_index( exo_index( arg_schema , Element ) )
{
}

FileSchema::~FileSchema() {}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

const FilePart * internal_declare_part(
  Schema            & arg_schema ,
  const std::string & arg_name ,
  int                 arg_id ,
  EntityType          arg_type ,
  ElementType         arg_element_type ,
  int                 arg_element_dim ,
  int                 arg_number_nodes ,
  int                 arg_number_attr ,
  const Field<double,1> & arg_attributes )
{
  static const char method[] = "phdmesh::exodus::FileSchema::declare_part" ;

  Part & part = arg_schema.declare_part( arg_name );

  CSet::Span<FilePart> attr = part.attribute<FilePart>();

  const FilePart * file_part = NULL ;

  if ( attr.empty() ) {

    const unsigned dim = arg_element_dim ;

    unsigned n = arg_number_attr ;

    if ( ! n ) { // Default number of attributes
      switch( arg_element_type ) {
      case CIRCLE :
      case SPHERE :
      case TRUSS :
      case SHELL : n = 1 ; break ;
      case BEAM :  n = dim == 2 ? 3 : ( dim == 3 ? 7 : 0 ); break ;
      default:     n = 0 ; break ;
      }
    }

    if ( n ) {
      arg_schema.declare_field_dimension(
        const_cast<Field<double,1>&>( arg_attributes ) , part , n );
    }

    file_part = new FilePart( part, arg_id, arg_type,
                              arg_element_type, arg_number_nodes , n );

    arg_schema.declare_part_attribute<FilePart>( part , file_part , true );
  }
  else if ( attr.size() == 1 ) {
    file_part = & *attr ;
  }

  if ( file_part == NULL ||
       & file_part->m_part         != & part ||
         file_part->m_identifier   != arg_id ||
         file_part->m_entity_type  != arg_type ||
         file_part->m_element_type != arg_element_type ||
         file_part->m_number_nodes != arg_number_nodes ||
         file_part->m_number_attr  != arg_number_attr ) {
    std::ostringstream msg ;
    msg << method << "FAILED redeclaration (count = " ;
    msg << attr.size();
    msg << ")" ;
    throw std::runtime_error( msg.str() );
  }

  return file_part ;
}

}

Part & FileSchema::declare_part(
  const std::string & arg_name ,
  int                 arg_id ,
  ElementType         arg_element_type ,
  unsigned            arg_number_nodes ,
  unsigned            arg_num_attributes )
{
  const FilePart * const fp =
    internal_declare_part( m_schema , arg_name , arg_id ,
                           Element , arg_element_type , m_dimension ,
                           arg_number_nodes , arg_num_attributes ,
                           m_field_elem_attr );

  m_parts[ Element ].push_back( fp );

  return fp->m_part ;
}

Part & FileSchema::declare_part(
  const std::string & arg_name ,
  int                 arg_id ,
  EntityType          arg_type )
{
  const FilePart * const fp =
     internal_declare_part( m_schema , arg_name , arg_id ,
                            arg_type , UNDEFINED , m_dimension , 0 , 0 ,
                            m_field_elem_attr );

  m_parts[ arg_type ].push_back( fp );

  return fp->m_part ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

struct BlockEntityProc {
  entity_id_type entity_id ;
  unsigned       block_id ;
  unsigned       processor ;
};

struct LessIndices {

  bool operator()( const std::pair<int,const Entity*> & L ,
                   const std::pair<int,const Entity*> & R ) const
  {
    return L.first != R.first ? L.first < R.first :
           L.second->identifier() < R.second->identifier() ;
  }

  bool operator ()( const BlockEntityProc & L ,
                    const BlockEntityProc & R ) const
  {
     return L.block_id  != R.block_id ?  L.block_id  < R.block_id : (
            L.entity_id != R.entity_id ? L.entity_id < R.entity_id :
                                         L.processor < R.processor );
  }
};

void assign_contiguous_indices(
  Mesh & M ,
  const std::vector< std::pair<int,const Entity*> > & entities ,
  const Field<int,1> & f )
{
  static const char method[] = "phdmesh::exodus::assign_contiguous_indices" ;

  const unsigned u_zero = 0 ;

  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();
  const ParallelMachine p_comm = M.parallel();

  std::vector< std::pair<int,const Entity*> >::const_iterator i ;

  unsigned unique_count = 0 ;

  // Get the total count, split evenly among processors, and then map
 
  entity_id_type end_id = 1 ;
  {
    const EntitySet & es = M.entities( f.entity_type() );
    if ( ! es.empty() ) {
      end_id += (-- es.end() )->identifier();
    }
  }

  all_reduce( p_comm , ReduceMax<1>( & end_id ) );

  const double p_map = ((double) p_size) / ((double) end_id) ;

  CommAll all( p_comm );

  for ( i = entities.begin() ; i != entities.end() ; ++i ) {
    const unsigned p = (unsigned)( p_map * i->second->identifier() );
    all.send_buffer( p ).skip<entity_id_type>(2);
  }

  all.allocate_buffers( p_size / 4 );

  for ( i = entities.begin() ; i != entities.end() ; ++i ) {
    entity_id_type data[2] ;
    data[0] = i->first ;
    data[1] = i->second->identifier();

    const unsigned p = (unsigned)( p_map * data[1] );
    all.send_buffer( p ).pack<entity_id_type>( data , 2);
  }

  all.communicate();

  std::vector<BlockEntityProc> received ;

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {
      entity_id_type data[2] ;
      buf.unpack<entity_id_type>( data , 2 );
      BlockEntityProc tmp ;
      tmp.block_id  = data[0] ;
      tmp.entity_id = data[1] ;
      tmp.processor = p ;
      received.push_back( tmp );
    }
  }

  std::sort( received.begin() , received.end() , LessIndices() );

  unique_count = 0 ;

  {
    entity_id_type id = 0 ;

    for ( std::vector<BlockEntityProc>::iterator
          j = received.begin() ; j != received.end() ; ++j ) {
      if ( id != j->entity_id ) { id = j->entity_id ; ++unique_count ; }
    }
  }

  std::vector<unsigned> global_count( p_size , u_zero );

  allgather( p_comm , & unique_count , & global_count[0] , 1 );

  entity_id_type count = 0 ;

  for ( unsigned p = 0 ; p < p_rank ; ++p ) { count += global_count[p] ; }

  entity_id_type global_index = count ;

  for ( unsigned p = p_rank ; p < p_size ; ++p ) { count += global_count[p] ; }

  all.swap_send_recv();
  all.reset_buffers();

  {
    entity_id_type id = 0 ;

    for ( std::vector<BlockEntityProc>::iterator
          j = received.begin() ; j != received.end() ; ++j ) {

      if ( id != j->entity_id ) { id = j->entity_id ; ++global_index ; }

      entity_id_type data[2] ;
      data[0] = global_index ;
      data[1] = id ;
      all.send_buffer( j->processor ).pack<entity_id_type>( data , 2 );
    }
  }

  all.communicate();

  entity_id_type local_index = 0 ;

  for ( i = entities.begin() ; i != entities.end() ; ++i ) {
    const entity_id_type id = i->second->identifier();

    const unsigned p = (unsigned)( p_map * id );
    CommBuffer & buf = all.recv_buffer( p );
    entity_id_type data[2] ;
    buf.unpack<entity_id_type>( data , 2 );

    int * const entity_index = i->second->data(f);

    entity_index[0] = data[0] ; // global index
    entity_index[1] = ++local_index ; // local index

    if ( id != data[1] ) {
      std::ostringstream msg ;
      msg << method << " FAILED" ;
      throw std::runtime_error(msg.str());
    }
  }
}

void assign_contiguous_indices( Mesh & M , const Field<int,1> & f )
{
  const EntitySet & es = M.entities( f.entity_type() );

  std::vector< std::pair<int,const Entity*> > entities ;

  entities.reserve( es.size() );

  std::pair<int,const Entity*> tmp ; tmp.first = 0 ;

  for ( EntitySet::const_iterator j = es.begin() ; j != es.end() ; ++j ) {
    tmp.second = & *j ;
    entities.push_back( tmp );
  }

  assign_contiguous_indices( M , entities , f );
}

}


void FileSchema::assign_indices( Mesh & arg_mesh ) const
{
  static const char method[] = "phdmesh::FileSchema::assign_indices" ;

  Mesh & M = arg_mesh ;
  const Schema & S = m_schema ;

  S.assert_same_schema( method , M.schema() );

  //--------------------------------------------------------------------

  assign_contiguous_indices( M , m_field_node_index );
  assign_contiguous_indices( M , m_field_edge_index );
  assign_contiguous_indices( M , m_field_face_index );

  //--------------------------------------------------------------------
  // Assign element indices by block and then element id order.
  // Count the number of owned elements per block.

  std::vector< std::pair<int,const Entity*> > entities ;

  entities.reserve( M.entities( Element ).size() );

  const std::vector<const FilePart*> & elem_parts = m_parts[ Element ];

  const unsigned num_elem_blk = elem_parts.size();

  const KernelSet & ks = M.kernels( Element );

  for ( KernelSet::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {

    const Kernel & k = *ik ;

    bool found = false ;

    for ( unsigned i = 0 ; i < num_elem_blk ; ++i ) {

      if ( k.has_superset( elem_parts[i]->m_part ) ) {

        if ( found ) {
          std::ostringstream msg ;
          msg << method
              << " FAILED: element contained in multiple blocks" ;
          throw std::runtime_error( msg.str() );
        }

        found = true ;

        std::pair<int,const Entity*> tmp ;
        tmp.first = elem_parts[i]->m_identifier ;

        for ( unsigned j = 0 ; j < k.size() ; ++j ) {
          tmp.second = k[j] ;
          entities.push_back( tmp );
        }
      }
    }
  }

  std::sort( entities.begin() , entities.end() , LessIndices() );

  // Assign element indices in element block order and then element id order

  assign_contiguous_indices( M , entities , m_field_elem_index );
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#if defined( PHDMESH_HAS_SNL_EXODUSII )

#include <exodusII.h>
#include <ne_nemesisI.h>

#if defined( PHDMESH_HAS_MPI )

namespace {

void broadcast( phdmesh::ParallelMachine c , int p_root , int * value , int n )
{ MPI_Bcast( value , n , MPI_INT , p_root , c ); }

void broadcast( phdmesh::ParallelMachine c , int p_root , char * value , int n )
{ MPI_Bcast( value , n , MPI_BYTE , p_root , c ); }

void scatter( phdmesh::ParallelMachine c ,
              const unsigned p_root ,
              const unsigned * const send_size ,
              const unsigned recv_size ,
              void * const data )
{
  const unsigned p_rank = phdmesh::parallel_machine_rank( c );
  const unsigned p_size = phdmesh::parallel_machine_size( c );

  std::vector<int> count ;
  std::vector<int> displ ;

  if ( p_rank == p_root ) {
    count.resize( p_size );
    displ.resize( p_size );
    unsigned n = 0 ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      count[p] = p_root != p ? send_size[p] : 0 ;
      displ[p] = n ;
      n += send_size[p] ;
    }
  }

  int * const p_count = p_rank == p_root ? & count[0] : NULL ;
  int * const p_displ = p_rank == p_root ? & displ[0] : NULL ;
  int p_recv = p_rank == p_root ? 0 : recv_size ;

  MPI_Scatterv( data , p_count , p_displ , MPI_BYTE ,
                data , p_recv ,            MPI_BYTE ,
                p_root , c );

  // Now the local copy, if any

  if ( p_rank == p_root && send_size[ p_root ] && displ[ p_root ] ) {
    unsigned char * const dst = reinterpret_cast<unsigned char *>( data );
    unsigned char * const src = dst + displ[ p_root ] ;
    memmove( dst , src , send_size[ p_root ] );
  }
}

}

#else

namespace {

void broadcast( phdmesh::ParallelMachine , int , char * , int ) {}
void broadcast( phdmesh::ParallelMachine , int , int * , int ) {}

void scatter( phdmesh::ParallelMachine ,
              const unsigned ,
              const unsigned * const ,
              const unsigned ,
              void * const ) {}

}

#endif


namespace phdmesh {
namespace exodus {

enum { Maximum_Name_Length = MAX_STR_LENGTH };

namespace {

const char name_circle[]   = "CIRCLE" ;
const char name_sphere[]   = "SPHERE" ;
const char name_truss[]    = "TRUSS" ;
const char name_beam[]     = "BEAM" ;
const char name_shell[]    = "SHELL" ;
const char name_quad[]     = "QUAD" ;
const char name_triangle[] = "TRIANGLE" ;
const char name_pyramid[]  = "PYRAMID" ;
const char name_tetra[]    = "TETRA" ;
const char name_wedge[]    = "WEDGE" ;
const char name_hex[]      = "HEX" ;
const char name_undefined[] = "UNDEFIEND" ;

const char * element_labels[] = {
  name_circle , name_sphere , name_truss , name_beam , name_shell ,
  name_quad , name_triangle ,
  name_pyramid , name_tetra , name_wedge , name_hex ,
  name_undefined
};

ElementType element_type( const char * label )
{
  ElementType type = UNDEFINED ;

  for ( unsigned i = 0 ; type == UNDEFINED && i < UNDEFINED ; ++i ) {
    const unsigned n = strlen( element_labels[i] );
    if ( ! strncmp( element_labels[i] , label , n ) ) {
      type = ElementType(i);
    }
  }
  return type ;
}

}

//----------------------------------------------------------------------

FileOutput::~FileOutput()
{
  const Mesh       & M  = m_mesh ;
  const FileSchema & FS = m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = FS.m_io_rank ;

  if ( p_rank == p_write ) { ex_close( m_exo_id ); }
}

//----------------------------------------------------------------------

namespace {

template<class Op>
EntitySet::const_iterator
  chunk( EntitySet::const_iterator i ,
         const EntitySet::const_iterator i_end ,
         Part & part ,
         const Field<int,1> & field_index ,
         const int index_begin ,
         const int index_end ,
         Op & op )
{
  const int maximum = op.max();

  std::vector<const Entity *> send_entities ;

  int index_begin_chunk = index_begin ;

  while ( index_begin_chunk < index_end ) {

    const int index_end_chunk =
      index_end < index_begin_chunk + maximum ?
      index_end : index_begin_chunk + maximum ;

    send_entities.clear();

    for ( bool flag = true ; flag && i != i_end ; ) {
      const Entity * const e = & *i ;
      if ( e->kernel().has_superset( part ) ) {
        flag = e->data(field_index)[0] < index_begin_chunk ;
      }
      if ( flag ) { ++i ; }
    }

    for ( bool flag = true ; flag && i != i_end ; ) {
      const Entity * const e = & *i ;
      if ( e->kernel().has_superset( part ) ) {
        if ( flag = e->data(field_index)[0] < index_end_chunk ) {
          send_entities.push_back( e );
        }
      }
      if ( flag ) { ++i ; }
    }

    op( index_begin_chunk , index_end_chunk , send_entities );

    index_begin_chunk = index_end_chunk ;
  }

  return i ;
}

//----------------------------------------------------------------------

struct WriteNodeIndexCoord {
  // Two buffers: communication and file data
  enum { size_per_node = 2 * ( sizeof(int) + 3 * sizeof(double) ) };

  FileOutput & output ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_error ;
  const char * exo_func ;
  std::vector<int>    identifiers ;
  std::vector<double> coord_x_y_z ;

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *> & );

  WriteNodeIndexCoord( FileOutput & arg_output , unsigned );
};

WriteNodeIndexCoord::WriteNodeIndexCoord(
  FileOutput & arg_output , unsigned max_buffer )
  : output( arg_output ),
    maximum( max_buffer / size_per_node ),
    writer( false ),
    exo_error( 0 ),
    exo_func( NULL ),
    identifiers(),
    coord_x_y_z()
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  writer = p_rank == p_write ;
}

void WriteNodeIndexCoord::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *> & send )
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_size =
     ( sizeof(int) * 2 + sizeof(double) * 3 ) * send.size();

  CommGather gather( p_comm , p_write , send_size );

  {
    CommBuffer & buf = gather.send_buffer();

    for ( std::vector<const Entity *>::const_iterator
          i = send.begin() ; i != send.end() ; ++i ) {
      const Entity & node = **i ;
      int ids[2] ;
      ids[0] = * node.data( SF.m_field_node_index );
      ids[1] =   node.identifier();
      const double * const xyz = node.data( SF.m_field_node_coord );

      buf.pack<int>( ids , 2 );
      buf.pack<double>( xyz , 3 );
    }
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id = output.exo_id() ;
    const int number = index_end - index_begin ;
    const unsigned n = number ;

    if ( identifiers.size() < n ) { identifiers.resize( n ); }
    if ( coord_x_y_z.size() < n ) { coord_x_y_z.resize( n * 3 ); }

    double * const x = & coord_x_y_z[0] ;
    double * const y = x + n ;
    double * const z = y + n ;

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = gather.recv_buffer(p);
      while ( buf.remaining() ) {
        int ids[2] ;    buf.unpack<int>( ids , 2 );
        double xyz[3] ; buf.unpack<double>( xyz , 3 );

        const unsigned offset = ids[0] - index_begin ;
        identifiers[ offset ] = ids[1] ;
        x[ offset ] = xyz[0] ;
        y[ offset ] = xyz[1] ;
        z[ offset ] = xyz[2] ;
      }
    }

    if ( ! exo_error ) {
      exo_func = "ne_put_n_node_num_map" ;
      exo_error =
        ne_put_n_node_num_map( exo_id, index_begin, number, & identifiers[0] );
    }

    if ( ! exo_error ) {
      exo_func = "ne_put_n_coord" ;
      exo_error = ne_put_n_coord(exo_id,index_begin,number,x,y,z);
    }
  }
}

//----------------------------------------------------------------------

struct WriteElemIndexConnect {
  FileOutput & output ;
  const FilePart & part ;
  int          index_part ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_error ;
  const char * exo_func ;

  std::vector<int> send_data ;
  std::vector<int> identifiers ;
  std::vector<int> connectivity ;

  WriteElemIndexConnect( FileOutput & ,
                         const FilePart & ,
                         int ,
                         unsigned max_buffer );

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *> & );
};

WriteElemIndexConnect::WriteElemIndexConnect(
  FileOutput & arg_output ,
  const FilePart   & arg_part ,
  int          arg_index_part ,
  unsigned max_buffer )
  : output( arg_output ), part( arg_part ),
    index_part( arg_index_part ), maximum( 0 ),
    writer( false ), exo_error( 0 ), exo_func( NULL ),
    send_data(), identifiers(), connectivity()
{
  const int i_zero = 0 ;
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_data_num = part.m_number_nodes + 2 ;

  // Two buffers: communication and file data

  maximum = max_buffer / ( 2 * send_data_num * sizeof(int) );

  send_data.resize( send_data_num , i_zero );

  writer = p_rank == p_write ;
}

void WriteElemIndexConnect::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *> & send )
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned block_id           = part.m_identifier ;
  const unsigned num_nodes_per_elem = part.m_number_nodes ;
  const unsigned send_data_num      = num_nodes_per_elem + 2 ;

  const unsigned send_size = sizeof(int) * send_data_num * send.size();

  CommGather gather( p_comm , p_write , send_size );

  {
    CommBuffer & buf = gather.send_buffer();

    for ( std::vector<const Entity *>::const_iterator
          i = send.begin() ; i != send.end() ; ++i ) {
      const Entity & elem = **i ;

      send_data[0] = * elem.data( SF.m_field_elem_index );
      send_data[1] =   elem.identifier();

      ConnectSpan con = elem.connections();

      for ( unsigned j = 0 ; j < num_nodes_per_elem ; ++j , ++con ) {
        Entity & node = * con->entity();
        send_data[j+2] = * node.data( SF.m_field_node_index );
      }

      buf.pack<int>( & send_data[0] , send_data_num );
    }
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id = output.exo_id() ;
    const int number = index_end - index_begin ;
    const unsigned n = number ;

    if ( identifiers.size() < n ) { identifiers.resize( n ); }
    if ( connectivity.size() < n * part.m_number_nodes ) {
      connectivity.resize( n * part.m_number_nodes );
    }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {

      CommBuffer & buf = gather.recv_buffer(p);

      while ( buf.remaining() ) {

        buf.unpack<int>( & send_data[0] , send_data_num );

        int index = send_data[0] ;

        unsigned offset = index - index_begin ;

        identifiers[offset] = send_data[1] ;

        int * const node_id = & connectivity[ offset * num_nodes_per_elem ];

        for ( unsigned j = 0 ; j < num_nodes_per_elem ; ++j ) {
          node_id[j] = send_data[j+2] ;
        }
      }
    }

    if ( ! exo_error ) {
      exo_func = "ne_put_n_elem_num_map" ;
      exo_error =
        ne_put_n_elem_num_map(exo_id,index_begin,number, & identifiers[0] );
    }

    if ( ! exo_error ) {
      const int index_begin_part = 1 + index_begin - index_part ;

      exo_func = "ne_put_n_elem_conn" ;
      exo_error =
        ne_put_n_elem_conn( exo_id , block_id ,
                            index_begin_part , number , & connectivity[0] );
    }
  }
}

unsigned count( const KernelSet & ks , const PartSet & ps )
{
  unsigned n = 0 ;

  KernelSet::const_iterator i ;
  for ( i = ks.begin() ; i != ks.end() ; ++i ) {
    if ( i->has_superset( ps ) ) {
      n += i->size();
    }
  }
  return n ;
}

std::string variable_name( const Part & p , 
                           const Field<void,0> & f ,
                           const unsigned k )
{
  std::string name ;

  const FieldDimension & dim = f.dimension( p );

  if ( 1 < dim.length() ) {
    unsigned dimension[ MaximumFieldDimension ];
    unsigned indices[   MaximumFieldDimension ];

    dim.dimension( dimension );
    dim.indices( k , indices );

    const FieldName * des = io_get_field_name( f );

    if ( des == NULL ) { des = & io_array_name(); }

    name = des->encode( f.name() , dim.number_of_dimensions() ,
                        dimension , indices );
  }
  else {
    name = f.name();
  }
  return name ;
}


void variable_add( const Field<void,0>                    & field ,
                   const std::vector<const FilePart *>    & parts ,
                   std::vector< std::string >             & var_names ,
                   std::vector< FieldIO > & spec )
{
  const unsigned size_begin = var_names.size();

  for ( unsigned j = 0 ; j < parts.size() ; ++j ) {
    Part & p = parts[j]->m_part ;

    const FieldDimension & d = field.dimension( p );

    for ( unsigned k = 0 ; k < d.length() ; ++k ) {
      const std::string name( variable_name( p , field , k ) );
      var_names.push_back( name );
    }
  }

  std::vector<std::string>::iterator i_beg = var_names.begin();

  std::advance( i_beg , size_begin );

  std::sort( i_beg , var_names.end() );

  var_names.erase( std::unique( i_beg , var_names.end() ) ,
                                        var_names.end() );
  
  FieldIO tmp ;

  tmp.m_field     = & field ;
  tmp.m_part      = NULL ;
  tmp.m_offset    = 0 ;
  tmp.m_var_index = 0 ;

  for ( unsigned j = 0 ; j < parts.size() ; ++j ) {
    tmp.m_part = parts[j] ;

    Part & p = parts[j]->m_part ;

    const FieldDimension & d = field.dimension( p );

    for ( unsigned k = 0 ; k < d.length() ; ++k ) {

      const std::string name( variable_name( p , field , k ) );

      const std::vector<std::string>::iterator i =
        std::lower_bound( i_beg , var_names.end() , name );

      tmp.m_offset = k ;
      tmp.m_var_index = 1 + std::distance( var_names.begin() , i );

      spec.push_back( tmp );
    }
  }
}

}

//----------------------------------------------------------------------

FileOutput::FileOutput(
  const FileSchema  & arg_schema ,
  const Mesh        & arg_mesh ,
  const std::string & arg_file_path ,
  const std::string & arg_title ,
  const bool          arg_storage_double ,
  const std::vector<const Field<void,0> * > & arg_fields ,
  const int * const arg_processor )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
{
  static const char method[] = "phdmesh::exodus::FileOutput::FileOutput" ;

  const int i_zero = 0 ;

  const Mesh       & M  = arg_mesh ;
  const Schema     & SM = M.schema();
  const FileSchema & FS = arg_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = FS.m_io_rank ;

  const bool writer = p_write == p_rank ;

  const Part & universal_part = SM.universal_part();
        Part & owns_part      = SM.owns_part();

  const std::vector<const FilePart *> & node_parts = FS.parts( Node );
  const std::vector<const FilePart *> & edge_parts = FS.parts( Edge );
  const std::vector<const FilePart *> & face_parts = FS.parts( Face );
  const std::vector<const FilePart *> & elem_parts = FS.parts( Element );

  //--------------------------------------------------------------------

  std::vector< std::string > name_node_var ;
  std::vector< std::string > name_elem_var ;
  std::vector<int> exist_elem_var ; 

  if ( arg_processor ) {
    const std::string name_processor("processor");

    FieldIO tmp ;
    tmp.m_field = NULL ;
    tmp.m_part  = NULL ;
    tmp.m_offset = 0 ;
    tmp.m_var_index = 1 ;

    if ( arg_processor[ Node ] ) {
      name_node_var.push_back( name_processor );
      m_field_node_universal.push_back( tmp );
    }

    if ( arg_processor[ Element ] ) {
      name_elem_var.push_back( name_processor );
      m_field_elem.push_back( tmp );
    }
  }

  for ( std::vector< const Field<void,0> * >::const_iterator
        i = arg_fields.begin() ; i != arg_fields.end() ; ++i ) {

    const Field<void,0>  & f = **i ;
    const FieldDimension & d = f.dimension( universal_part );

    switch( f.entity_type() ) {
    case Node :
      if ( d.length() ) { // Universal nodal variable
        FieldIO tmp ;

        tmp.m_field     = & f ;
        tmp.m_part      = NULL ;
        tmp.m_offset    = 0 ;
        tmp.m_var_index = 0 ;

        for ( unsigned k = 0 ; k < d.length() ; ++k ) {
          name_node_var.push_back( variable_name(universal_part, f, k) );
          tmp.m_offset    = k ;
          tmp.m_var_index = name_node_var.size();
          m_field_node_universal.push_back( tmp );
        }
      }
      else { // Node set variable
        // variable_add( f , node_parts , name_nodeset_var , m_fields );
      }
      break ;

    case Element :
      variable_add( f , elem_parts , name_elem_var , m_field_elem );
      break ;

    default :
      break ;
    }
  }

  exist_elem_var.resize( elem_parts.size() * name_elem_var.size() , i_zero);

  for ( unsigned j = 0 ; j < elem_parts.size() ; ++j ) {

    int * const exist = & exist_elem_var[ j * name_elem_var.size() ];

    for ( unsigned k = 0 ; k < m_field_elem.size() ; ++k ) {
      if ( m_field_elem[k].m_part == NULL ||
           m_field_elem[k].m_part == elem_parts[k] ) { 
        const int offset = m_field_elem[k].m_var_index - 1 ;
        exist[ offset ] = 1 ;
      }
    }
  }

  //--------------------------------------------------------------------
  // Sizes and counts

  const int num_dim       = FS.m_dimension ;
  const int num_node_sets = node_parts.size();
  const int num_edge_sets = edge_parts.size();
  const int num_face_sets = face_parts.size();
  const int num_side_sets = num_edge_sets + num_face_sets ;
  const int num_elem_blk  = elem_parts.size();

  const unsigned num_parts = num_node_sets + num_side_sets + num_elem_blk ;

  m_global_counts.resize( 4 + num_parts , i_zero );

  std::vector<int> part_count_local(  4 + num_parts , i_zero );

  {
    PartSet ps ;
    insert( ps , owns_part );

    part_count_local[ 0 ] = count( M.kernels( Node ) , ps );
    part_count_local[ 1 ] = count( M.kernels( Edge ) , ps );
    part_count_local[ 2 ] = count( M.kernels( Face ) , ps );
    part_count_local[ 3 ] = count( M.kernels( Element ) , ps );
  }

  unsigned n = 4 ;
  for ( int i = 0 ; i < num_node_sets ; ++i , ++n ) {
    Part & part = node_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Node ) , ps );
  }

  for ( int i = 0 ; i < num_edge_sets ; ++i , ++n ) {
    Part & part = edge_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Edge ) , ps );
  }

  for ( int i = 0 ; i < num_face_sets ; ++i , ++n ) {
    Part & part = face_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Face ) , ps );
  }

  for ( int i = 0 ; i < num_elem_blk ; ++i , ++n ) {
    Part & part = elem_parts[i]->m_part ;
    PartSet ps ;
    insert( ps , part );
    insert( ps , owns_part );
    part_count_local[ n ] = count( M.kernels( Element ) , ps );
  }

  all_reduce_sum( p_comm , & part_count_local[0] , & m_global_counts[0] , n );

  const int num_nodes_global = m_global_counts[0] ;
  const int num_elems_global = m_global_counts[3] ;

  const int * const global_num_elem_in_block =
    & m_global_counts[ 4 + num_node_sets + num_edge_sets + num_face_sets ];

  //--------------------------------------------------------------------

  int exo_error = 0 ;
  const char * exo_func = NULL ;

  if ( writer ) {

    // Create file:
    {
      int comp_ws = sizeof(double);
      int io_ws   = arg_storage_double ? sizeof(double) : sizeof(float) ;

      exo_func = "ex_create" ;
      m_exo_id = ex_create( arg_file_path.c_str() ,
                            EX_CLOBBER , & comp_ws , & io_ws );

      if ( m_exo_id < 0 ) { exo_error = -1 ; }
    }

    // Put title and sizes:
    if ( ! exo_error ) {
      exo_func = "ex_put_init" ;
      exo_error = ex_put_init( m_exo_id , arg_title.c_str() ,
                               num_dim , num_nodes_global , num_elems_global ,
                               num_elem_blk , num_node_sets , num_side_sets );
    }

    // Put the model nodal coordinate names:
    if ( ! exo_error ) {
      char coord_name_x[ Maximum_Name_Length ] ;
      char coord_name_y[ Maximum_Name_Length ] ;
      char coord_name_z[ Maximum_Name_Length ] ;

      char * coord_names[3] ;

      coord_names[0] = coord_name_x ;
      coord_names[1] = coord_name_y ;
      coord_names[2] = coord_name_z ;

      strcpy( coord_name_x , FS.m_field_node_coord.name().c_str() );
      strcpy( coord_name_y , FS.m_field_node_coord.name().c_str() );
      strcpy( coord_name_z , FS.m_field_node_coord.name().c_str() );
      strcat( coord_name_x , "_x" );
      strcat( coord_name_y , "_y" );
      strcat( coord_name_z , "_z" );

      exo_func = "ex_put_coord_names" ;
      exo_error = ex_put_coord_names( m_exo_id , coord_names );
    }

    // Put element blocks' description:
    if ( ! exo_error ) {
      char * const char_ptr_null = NULL ;

      std::vector<char> tmp( num_elem_blk * Maximum_Name_Length );

      std::vector<int>    elem_blk_id( num_elem_blk , 0 );
      std::vector<char *> elem_blk_name( num_elem_blk , char_ptr_null );
      std::vector<char *> elem_type( num_elem_blk , char_ptr_null );
      std::vector<int>    num_nodes_per_elem( num_elem_blk , 0 );
      std::vector<int>    num_attr( num_elem_blk , 0 );

      for ( int i = 0 ; i < num_elem_blk ; ++i ) {
        const FilePart   & fp = * elem_parts[i] ;
        elem_blk_id[i]        = fp.m_identifier ;
        num_nodes_per_elem[i] = fp.m_number_nodes ;
        num_attr[i]           = fp.m_number_attr ;
        elem_type[i]          = & tmp[ i * Maximum_Name_Length ] ;
        elem_blk_name[i]      = const_cast<char*>( fp.m_part.name().c_str() );

        strcpy( elem_type[i] , element_labels[ fp.m_element_type ] );
      }

      exo_func = "ex_put_concat_elem_block" ;
      exo_error = ex_put_concat_elem_block( m_exo_id ,
                                            & elem_blk_id[0] ,
                                            & elem_type[0] ,
                                              global_num_elem_in_block ,
                                            & num_nodes_per_elem[0] ,
                                            & num_attr[0] ,
                                            1 /* have identifier maps */ );

      if ( ! exo_error ) {
        exo_func = "ex_put_names(element)" ;
        exo_error = ex_put_names( m_exo_id, EX_ELEM_BLOCK, & elem_blk_name[0]);
      }
    }

    // Put results variable description:

    if ( ! exo_error ) {

      int num_var_global = 0 ;
      int num_var_node = name_node_var.size();
      int num_var_elem = name_elem_var.size();
      int num_var_side_set = 0 ;
      int num_var_node_set = 0 ;

      exo_func = "ex_put_all_var_param" ;
      exo_error = ex_put_all_var_param( m_exo_id ,
                                        num_var_global ,
                                        num_var_node ,
                                        num_var_elem ,
                                        & exist_elem_var[0] ,
                                        num_var_node_set ,
                                        NULL ,
                                        num_var_side_set ,
                                        NULL );

      if ( ! exo_error && num_var_node ) {
        std::vector<char *> ptr_var_names( name_node_var.size() );

        for ( unsigned i = 0 ; i < name_node_var.size() ; ++i ) {
          ptr_var_names[i] = const_cast<char*>( name_node_var[i].c_str() );
        }

        exo_func = "ex_put_var_names(node)" ;
        exo_error = ex_put_var_names( m_exo_id , "n" ,
                                      num_var_node ,
                                      & ptr_var_names[0] );
      }

      if ( ! exo_error && num_var_elem ) {
        std::vector<char *> ptr_var_names( name_elem_var.size() );

        for ( unsigned i = 0 ; i < name_elem_var.size() ; ++i ) {
          ptr_var_names[i] = const_cast<char*>( name_elem_var[i].c_str() );
        }

        exo_func = "ex_put_var_names(element)" ;
        exo_error = ex_put_var_names( m_exo_id , "e" ,
                                      num_var_elem ,
                                      & ptr_var_names[0] );
      }
    }
  }

  broadcast( p_comm , p_write , & exo_error , 1 );

  //--------------------------------------------------------------------
  // Write nodal identifiers and genesis coordinates.

  if ( ! exo_error ) {
    WriteNodeIndexCoord work( *this , m_max_buffer );

    const EntitySet & es = M.entities( Node );

    const int index = 1 ;
    const int index_end = index + num_nodes_global ;

    chunk( es.begin() , es.end() ,
           owns_part , FS.m_field_node_index ,
           index , index_end , work );

    exo_error = work.exo_error ;

    broadcast( p_comm , p_write , & exo_error , 1 );

    exo_func  = work.exo_func ;
  }

  //--------------------------------------------------------------------
  // Write element block identifiers and connectivity,
  // element indices are partitioned by element block.

  if ( ! exo_error ) {
    const EntitySet & es = M.entities( Element );

    EntitySet::const_iterator i = es.begin();

    int index = 1 ;

    for ( int j = 0 ; j < num_elem_blk ; ++j ) {

      const int index_end = index + global_num_elem_in_block[j] ;

      WriteElemIndexConnect work( *this, * elem_parts[j], index, m_max_buffer );

      i = chunk( i , es.end() , owns_part , FS.m_field_elem_index ,
                 index , index_end , work );

      index = index_end ;

      if ( ! exo_error ) {
        exo_error = work.exo_error ;
        exo_func  = work.exo_func ;
      }
    }

    broadcast( p_comm , p_write , & exo_error , 1 );
  }

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( writer ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void pack_entity_proc( const std::vector<const Entity *> & send ,
                       const Field<int,1>                & f_index ,
                       CommBuffer & buf )
{
  double item_buffer[2] ;

  for ( std::vector<const Entity *>::const_iterator
        i = send.begin() ; i != send.end() ; ++i ) {
    const Entity & e = **i ;

    item_buffer[0] = * e.data( f_index );
    item_buffer[1] =   e.owner_rank();

    buf.pack<double>( item_buffer , 2 );
  }
}

template<typename T>
void pack_entity_data( const std::vector<const Entity *> & send ,
                       const Field<int,1>                & f_index ,
                       const Field<T,0>                  & f_data ,
                       const unsigned                      offset ,
                       CommBuffer & buf )
{
  double item_buffer[2] ;

  for ( std::vector<const Entity *>::const_iterator
        i = send.begin() ; i != send.end() ; ++i ) {
    const Entity & e = **i ;

    item_buffer[0] = * e.data( f_index );
    item_buffer[1] =   e.data( f_data )[ offset ] ;

    buf.pack<double>( item_buffer , 2 );
  }
}

void pack_entity_data( const std::vector<const Entity *> & send ,
                       const Field<int,1>                & f_index ,
                       const Field<void,0>               & f_data ,
                       const unsigned                      offset ,
                       CommBuffer & buf )
{
  switch( f_data.numeric_type_ordinal() ) {
  case NumericEnum< signed char >::value :
    pack_entity_data( send, f_index, f_data.field<signed char,0>(), offset, buf );
    break ;
  case NumericEnum< unsigned char >::value :
    pack_entity_data( send, f_index, f_data.field<unsigned char,0>(), offset, buf );
    break ;
  case NumericEnum< signed short >::value :
    pack_entity_data( send, f_index, f_data.field<signed short,0>(), offset, buf );
    break ;
  case NumericEnum< unsigned short >::value :
    pack_entity_data( send, f_index, f_data.field<unsigned short,0>(), offset, buf );
    break ;
  case NumericEnum< signed int >::value :
    pack_entity_data( send, f_index, f_data.field<signed int,0>(), offset, buf );
    break ;
  case NumericEnum< unsigned int >::value :
    pack_entity_data( send, f_index, f_data.field<unsigned int,0>(), offset, buf );
    break ;
  case NumericEnum< signed long >::value :
    pack_entity_data( send, f_index, f_data.field<signed long,0>(), offset, buf );
    break ;
  case NumericEnum< unsigned long >::value :
    pack_entity_data( send, f_index, f_data.field<unsigned long,0>(), offset, buf );
    break ;
  case NumericEnum< float >::value :
    pack_entity_data( send, f_index, f_data.field<float,0>(), offset, buf );
    break ;
  case NumericEnum< double >::value :
    pack_entity_data( send, f_index, f_data.field<double,0>(), offset, buf );
    break ;
  default : break ;
  }
}

//----------------------------------------------------------------------

struct WriteNodeValues {

  FileOutput & output ;
  FieldIO * field ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_error ;
  const char * exo_func ;
  std::vector<double> recv_buffer ;

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *> & );

  WriteNodeValues( FileOutput & , unsigned );

  void set_field( FieldIO & f ) { field = & f ; }
};

WriteNodeValues::WriteNodeValues(
  FileOutput & arg_output , unsigned max_buffer )
  : output( arg_output ),
    field( NULL ),
    maximum( 0 ),
    writer( false ),
    exo_error( 0 ),
    exo_func( NULL ),
    recv_buffer()
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  writer = p_write == p_rank ;

  // Two buffers for field data and one for index

  maximum = max_buffer / ( 3 * sizeof(double) );
}

void WriteNodeValues::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *> & send )
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_size = 2 * sizeof(double) * send.size();

  double item_buffer[2] ;

  CommGather gather( p_comm , p_write , send_size );

  if ( field->m_field ) {
    pack_entity_data( send , SF.m_field_node_index ,
                      * field->m_field , field->m_offset ,
                       gather.send_buffer() );
  }
  else {
    pack_entity_proc( send , SF.m_field_node_index , gather.send_buffer() );
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id   = output.exo_id() ;
    const int exo_step = output.exo_step();
    const int number   = index_end - index_begin ;
    const unsigned n = number ;

    if ( recv_buffer.size() < n ) { recv_buffer.resize( n ); }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = gather.recv_buffer(p);
      while ( buf.remaining() ) {
        buf.unpack<double>( item_buffer , 2 );

        const unsigned index  = (unsigned) item_buffer[0] ;
        const unsigned offset = index - index_begin ;

        recv_buffer[ offset ] = item_buffer[1] ;
      }
    }

    exo_func = "ne_put_nodal_var_slab" ;
    exo_error = ne_put_nodal_var_slab( exo_id , exo_step, field->m_var_index ,
                                       index_begin , number ,
                                       & recv_buffer[0] );
  }
}

//----------------------------------------------------------------------

struct WriteElemValues {

  FileOutput & output ;
  FieldIO * field ;
  int          index_part ;
  unsigned     maximum ;
  bool         writer ;
  int          exo_ident ;
  int          exo_error ;
  const char * exo_func ;
  std::vector<double> recv_buffer ;

  unsigned max() const { return maximum ; }

  void operator()( const int , const int ,
                   const std::vector<const Entity *> & );

  WriteElemValues( FileOutput & , unsigned );

  void set_field( FieldIO & f , int index , int ident )
    { field = & f ; index_part = index ; exo_ident = ident ; }
};

WriteElemValues::WriteElemValues(
  FileOutput & arg_output , unsigned max_buffer )
  : output( arg_output ),
    field( NULL ),
    index_part( 0 ),
    maximum( 0 ),
    writer( false ),
    exo_ident( 0 ),
    exo_error( 0 ),
    exo_func( NULL ),
    recv_buffer()
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_write = SF.m_io_rank ;

  writer = p_write == p_rank ;

  // Two buffers for field data and one for index

  maximum = max_buffer / ( 3 * sizeof(double) );
}

void WriteElemValues::operator()(
  const int index_begin ,
  const int index_end ,
  const std::vector<const Entity *> & send )
{
  const Mesh       & M  = output.m_mesh ;
  const FileSchema & SF = output.m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_write = SF.m_io_rank ;

  const unsigned send_size = 2 * sizeof(double) * send.size();

  double item_buffer[2] ;

  CommGather gather( p_comm , p_write , send_size );

  if ( field->m_field ) {
    pack_entity_data( send , SF.m_field_elem_index ,
                      * field->m_field , field->m_offset ,
                       gather.send_buffer() );
  }
  else {
    pack_entity_proc( send , SF.m_field_elem_index , gather.send_buffer() );
  }

  gather.communicate();

  if ( writer ) {

    const int exo_id   = output.exo_id() ;
    const int exo_step = output.exo_step();
    const int number   = index_end - index_begin ;
    const unsigned n = number ;
    const int index_begin_part = 1 + index_begin - index_part ;

    if ( recv_buffer.size() < n ) { recv_buffer.resize( n ); }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = gather.recv_buffer(p);
      while ( buf.remaining() ) {
        buf.unpack<double>( item_buffer , 2 );

        const unsigned index  = (unsigned) item_buffer[0] ;
        const unsigned offset = index - index_begin ;

        recv_buffer[ offset ] = item_buffer[1] ;
      }
    }

    exo_func = "ne_put_elem_var_slab" ;
    exo_error = ne_put_elem_var_slab( exo_id , exo_step,
                                      field->m_var_index ,
                                      exo_ident ,
                                      index_begin_part , number ,
                                      & recv_buffer[0] );
  }
}

}

//----------------------------------------------------------------------

void FileOutput::write( double time_value )
{
  static const char method[] = "phdmesh::exodus::FieldIO::wrote" ;

  const Mesh       & M  = m_mesh ;
  const Schema     & SM = M.schema();
  const FileSchema & FS = m_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_rank = M.parallel_rank() ;
  const unsigned  p_write = FS.m_io_rank ;

  const std::vector<const FilePart *> & node_parts = FS.parts( Node );
  const std::vector<const FilePart *> & edge_parts = FS.parts( Edge );
  const std::vector<const FilePart *> & face_parts = FS.parts( Face );
  const std::vector<const FilePart *> & elem_parts = FS.parts( Element );

  const unsigned num_node_parts = node_parts.size();
  const unsigned num_edge_parts = edge_parts.size();
  const unsigned num_face_parts = face_parts.size();
  const unsigned num_elem_parts = elem_parts.size();

  const bool writer = p_write == p_rank ;

  Part & owns_part = SM.owns_part();

  int exo_error = 0 ;
  const char * exo_func = NULL ;

  ++m_counter ;

  if ( writer ) {
    exo_func = "ex_put_time" ;
    exo_error = ex_put_time( m_exo_id , m_counter , & time_value );
  }

  // exo_error is NOT parallel consistent

  {
    const int num_nodes_global = m_global_counts[ Node ] ;

    WriteNodeValues work( *this , m_max_buffer );

    for ( std::vector<FieldIO>::iterator
          i =  m_field_node_universal.begin() ;
          i != m_field_node_universal.end() ; ++i ) {

      work.set_field( *i );

      const EntitySet & es = M.entities( Node );

      const int index = 1 ;
      const int index_end = index + num_nodes_global ;

      chunk( es.begin() , es.end() ,
             owns_part , FS.m_field_node_index ,
             index , index_end , work );
    }

    exo_func  = work.exo_func ;
    exo_error = work.exo_error ;

    broadcast( p_comm , p_write , & exo_error , 1 );
  }

  if ( ! exo_error ) {
    // Iteration: Variables, element blocks, element chunking

    const EntitySet & es = M.entities( Element );

    const int * const global_num_elem_in_block =
      & m_global_counts[ 4 + num_node_parts + num_edge_parts + num_face_parts ];

    WriteElemValues work( *this , m_max_buffer );

    for ( std::vector<FieldIO>::iterator
          k =  m_field_elem.begin() ;
          k != m_field_elem.end() ; ++k ) {

      const FilePart * const part = k->m_part ;

      EntitySet::const_iterator i = es.begin();

      int index = 1 ;

      for ( unsigned j = 0 ; j < num_elem_parts ; ++j ) {

        const int index_end = index + global_num_elem_in_block[j] ;

        if ( part == NULL || part == elem_parts[j] ) {

          work.set_field( *k , index , elem_parts[j]->m_identifier );

          i = chunk( i , es.end() ,
                     owns_part , FS.m_field_elem_index ,
                     index , index_end , work );

          if ( ! exo_error ) {
            exo_error = work.exo_error ;
            exo_func  = work.exo_func ;
          }
        }
        else {
          // Skip elements in this block
          for ( bool flag = true ; flag && i != es.end() ; ) {
            if ( i->kernel().has_superset( owns_part ) ) {
              flag = i->data(FS.m_field_elem_index)[0] < index_end ;
            }
            if ( flag ) { ++i ; }
          }
        }
        index = index_end ;
      }
    }

    exo_func  = work.exo_func ;
    exo_error = work.exo_error ;
    broadcast( p_comm , p_write , & exo_error , 1 );
  }

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( p_rank == p_write ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Open and read file for meta-data

namespace {





enum { MAX_NUM_ATTR = 8 };

struct exo_elem_block_data {
  char name[ Maximum_Name_Length ];
  char type[ Maximum_Name_Length ];
  char attr[ MAX_NUM_ATTR ][ Maximum_Name_Length ];
  int  block_id ;
  int  num_nodes ;
  int  num_attr ;
  int  error ;
};

}

FileSchema::FileSchema( Schema & arg_schema ,
                        const Field<double,1> & arg_node_coordinates ,
                        const Field<double,1> & arg_elem_attributes ,
                        const std::string     & arg_file_path ,
                        ParallelMachine         arg_comm ,
                        const unsigned          arg_reader_rank )
  : m_schema( arg_schema ),
    m_io_rank( arg_reader_rank ),
    m_dimension( arg_node_coordinates.max_length() ),
    m_field_node_coord( arg_node_coordinates ),
    m_field_elem_attr(  arg_elem_attributes ),
    m_field_node_index( exo_index( arg_schema , Node ) ),
    m_field_edge_index( exo_index( arg_schema , Edge ) ),
    m_field_face_index( exo_index( arg_schema , Face ) ),
    m_field_elem_index( exo_index( arg_schema , Element ) )
{
  static const char method[] = "phdmesh::exodus::FileSchema::FileSchema" ;

  ParallelMachine p_comm = arg_comm ;
  const unsigned  p_rank = parallel_machine_rank( arg_comm );
  const unsigned  p_read = arg_reader_rank ;

  //--------------------------------------------------------------------

  int exo_data[ 32 ];

  int exo_id = 0 ;
  int exo_error = 0 ;
  const char * exo_func = NULL ;

  const bool reader = p_rank == p_read ;

  if ( reader ) { // Open file:
    int comp_ws = sizeof(double) ;
    int io_ws   = 0 ;
    float version = 0 ;

    exo_func = "ex_open" ;
    exo_id = ex_open( arg_file_path.c_str() , EX_READ ,
                        & comp_ws , & io_ws , & version );

    if ( exo_id < 0 ) { exo_error = -1 ; }
  }

  // Get sizes:

  if ( reader && ! exo_error ) {
    char title[ MAX_LINE_LENGTH ];
    exo_func = "ex_get_init" ;
    exo_error = ex_get_init( exo_id , title ,
                             exo_data + 1 ,
                             exo_data + 2 ,
                             exo_data + 3 ,
                             exo_data + 4 ,
                             exo_data + 5 ,
                             exo_data + 6 );
  }

  exo_data[0] = exo_error ;

  broadcast( p_comm , p_read , exo_data , 7 );

  exo_error = exo_data[0] ;

  const int num_dim          = exo_data[1] ;
  // const int num_nodes_global = exo_data[2] ;
  // const int num_elems_global = exo_data[3] ;
  const int num_elem_blk     = exo_data[4] ;
  // const int num_node_sets    = exo_data[5] ;
  // const int num_side_sets    = exo_data[6] ;

  if ( num_dim != (int) m_dimension ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED: incompatible spatial dimensions "
          << m_dimension << " != "
          << num_dim ;
    }
    throw std::runtime_error( msg.str() );
  }

  //--------------------------------------------------------------------
  // Get element blocks names and declare element parts

  if ( ! exo_error ) {

    std::vector<exo_elem_block_data> data( num_elem_blk );

    {
      char * begin = reinterpret_cast<char *>( & data[0] );
      char * end   = reinterpret_cast<char *>( & data[ num_elem_blk ] );
      const unsigned n = end - begin ;

      memset( begin , 0 , n );
    }

    if ( reader ) {
      const unsigned num_names =
        num_elem_blk < MAX_NUM_ATTR ? MAX_NUM_ATTR : num_elem_blk ;

      std::vector<char*> names( num_names );

      std::vector<int> ids( num_elem_blk );

      exo_func = "ex_get_elem_blk_ids" ;
      exo_error = ex_get_elem_blk_ids( exo_id , & ids[0] );

      if ( ! exo_error ) {

        for ( int i = 0 ; i < num_elem_blk ; ++i ) {
          data[i].block_id = ids[i] ;
          names[i] = data[i].name ;
        }

        exo_func = "ex_get_names" ;
        exo_error = ex_get_names( exo_id , EX_ELEM_BLOCK , & names[0] );

        if ( 0 < exo_error ) { // Names are not defined
          for ( int i = 0 ; i < num_elem_blk ; ++i ) {
            sprintf( data[i].name , "block_%d" , ids[i] );
          }
          exo_error = 0 ;
        }
      }

      for ( int i = 0 ; ! exo_error && i < num_elem_blk ; ++i ) {
        int num_elem_this_blk ; // discard for now

        exo_func = "ex_get_elem_block" ;
        exo_error = ex_get_elem_block( exo_id , ids[i] ,
                                       data[i].type ,
                                       & num_elem_this_blk ,
                                       & data[i].num_nodes ,
                                       & data[i].num_attr );

        if ( ! exo_error ) {
          for ( unsigned j = 0 ; j < MAX_NUM_ATTR ; ++j ) {
            names[i] = data[i].attr[j] ;
          }
          exo_func = "ex_get_elem_attr_names" ;
          exo_error = ex_get_elem_attr_names( exo_id , ids[i] , & names[0] );

          if ( 0 < exo_error ) {
            for ( unsigned j = 0 ; j < MAX_NUM_ATTR ; ++j ) {
              data[i].attr[j][0] = 0 ;
            }
            exo_error = 0 ;
          }
        }
      }

      data[0].error = exo_error ;
    }

    {
      char * begin = reinterpret_cast<char *>( & data[0] );
      char * end   = reinterpret_cast<char *>( & data[ num_elem_blk ] );
      const unsigned n = end - begin ;

      broadcast( p_comm , p_read , begin , n );

      exo_error = data[0].error ;
    }

    if ( ! exo_error ) {
      for ( int i = 0 ; i < num_elem_blk ; ++i ) {

        const FilePart * fp =
          internal_declare_part( m_schema ,
                                 std::string( data[i].name ) ,
                                 data[i].block_id ,
                                 Element ,
                                 element_type( data[i].type ) ,
                                 m_dimension ,
                                 data[i].num_nodes ,
                                 data[i].num_attr ,
                                 m_field_elem_attr );

        m_parts[ Element ].push_back( fp );
      }
    }
  }

  //--------------------------------------------------------------------

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }

  if ( reader ) { ex_close( exo_id ); }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// map index to processor:  p = ( p_size * i ) / n_total

unsigned index_processor_map( const unsigned p_size ,
                              const unsigned i ,
                              const unsigned n_total )
{
  return ( p_size * i ) / n_total ;
}

unsigned index_processor_count( const unsigned p ,
                                const unsigned p_size ,
                                const unsigned i ,
                                const unsigned n ,
                                const unsigned n_total )
{
  // The beginning and ending indices for the given processor
  // within the span of indices [i,i+n)
  unsigned i_beg = ( n_total * p ) / p_size ;
  unsigned i_end = ( n_total * ( p + 1 ) ) / p_size ;
  if ( i_beg < i ) { i_beg = i ; }
  if ( i + n < i_end ) { i_end = i + n ; }
  return i_beg < i_end ? ( i_end - i_beg ) : 0 ;
}

struct NodeData {
  double coord[3] ;
  int    ident ;
  int    index ;

  static unsigned size_of()
    {
      NodeData * tmp = NULL ;
      return reinterpret_cast<char*>(tmp+1) - reinterpret_cast<char*>(tmp);
    }
};

struct less_NodeData {
  bool operator()( const NodeData & lhs , const int rhs ) const
    { return lhs.index < rhs ; }
};

}

FileInput::~FileInput()
{
  const Mesh       & M  = m_mesh ;
  const FileSchema & FS = m_schema ;
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_read = FS.m_io_rank ;

  if ( p_rank == p_read ) { ex_close( m_exo_id ); }
}

FileInput::FileInput(
  const FileSchema  & arg_schema ,
        Mesh        & arg_mesh ,
  const std::string & arg_file_path ,
  const std::vector< const Field<void,0> * > & arg_fields )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
{
  static const char method[] = "phdmesh::exodus::FileInput::FileInput" ;

        Mesh       & M  = arg_mesh ;
  const Schema     & SM = M.schema();
  const FileSchema & FS = arg_schema ;
  ParallelMachine p_comm = M.parallel();
  const unsigned  p_size = M.parallel_size();
  const unsigned  p_rank = M.parallel_rank();
  const unsigned  p_read = FS.m_io_rank ;

  const bool reader = p_read == p_rank ;

  const Part & universal_part = SM.universal_part();
        Part & owns_part      = SM.owns_part();

  // const std::vector<const FilePart *> & node_parts = FS.parts( Node );
  // const std::vector<const FilePart *> & edge_parts = FS.parts( Edge );
  // const std::vector<const FilePart *> & face_parts = FS.parts( Face );
  const std::vector<const FilePart *> & elem_parts = FS.parts( Element );

  //--------------------------------------------------------------------

  std::vector< std::string > name_node_var ;
  std::vector< std::string > name_elem_var ;
  std::vector<int> exist_elem_var ; 

  for ( std::vector< const Field<void,0> * >::const_iterator
        i = arg_fields.begin() ; i != arg_fields.end() ; ++i ) {

    const Field<void,0>  & f = **i ;
    const FieldDimension & d = f.dimension( universal_part );

    switch( f.entity_type() ) {
    case Node :
      if ( d.length() ) { // Universal nodal variable
        FieldIO tmp ;

        tmp.m_field     = & f ;
        tmp.m_part      = NULL ;
        tmp.m_offset    = 0 ;
        tmp.m_var_index = 0 ;

        for ( unsigned k = 0 ; k < d.length() ; ++k ) {
          name_node_var.push_back( variable_name(universal_part, f, k) );
          tmp.m_offset    = k ;
          tmp.m_var_index = name_node_var.size();
          m_field_node_universal.push_back( tmp );
        }
      }
      else { // Node set variable
        // variable_add( f , node_parts , name_nodeset_var , m_fields );
      }
      break ;

    case Element :
      variable_add( f , elem_parts , name_elem_var , m_field_elem );
      break ;

    default :
      break ;
    }
  }

  //--------------------------------------------------------------------
  // Sizes and counts

  int exo_data[ 32 ];

  int exo_error = 0 ;
  const char * exo_func = NULL ;

  if ( reader ) { // Open file:
    int comp_ws = sizeof(double) ;
    int io_ws   = 0 ;
    float version = 0 ;

    exo_func = "ex_open" ;
    m_exo_id = ex_open( arg_file_path.c_str() , EX_READ ,
                        & comp_ws , & io_ws , & version );

    if ( m_exo_id < 0 ) { exo_error = -1 ; }
  }

  // Get sizes:

  if ( reader && ! exo_error ) {
    char title[ MAX_LINE_LENGTH ];
    exo_func = "ex_get_init" ;
    exo_error = ex_get_init( m_exo_id , title ,
                             exo_data + 1 ,
                             exo_data + 2 ,
                             exo_data + 3 ,
                             exo_data + 4 ,
                             exo_data + 5 ,
                             exo_data + 6 );

    if ( 0 == exo_error ) {
      int tmp ;

      exo_func = "ne_get_n_node_num_map" ;
      exo_error = ne_get_n_node_num_map( m_exo_id, 1, 1, & tmp );

      exo_data[7] = 0 == exo_error ;

      if ( 0 < exo_error ) { exo_error = 0 ; }
    }

    if ( 0 == exo_error ) {
      int tmp ;

      exo_func = "ne_get_n_elem_num_map" ;
      exo_error = ne_get_n_elem_num_map( m_exo_id, 1, 1, & tmp );

      exo_data[8] = 0 == exo_error ;

      if ( 0 < exo_error ) { exo_error = 0 ; }
    }
  }

  exo_data[0] = exo_error ;

  broadcast( p_comm , p_read , exo_data , 9 );

  exo_error = exo_data[0] ;

  const int num_dim          = exo_data[1] ;
  const int num_nodes_global = exo_data[2] ;
  const int num_elems_global = exo_data[3] ;
  const int num_elem_blk     = exo_data[4] ;
  // const int num_node_sets    = exo_data[5] ;
  // const int num_side_sets    = exo_data[6] ;
  const bool has_node_identifiers = exo_data[7] ;
  const bool has_elem_identifiers = exo_data[8] ;

  if ( num_dim != (int) FS.m_dimension ||
       num_elem_blk != (int) elem_parts.size() ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED: incompatible spatial dimension "
          << FS.m_dimension << " != "
          << num_dim
          << " OR element block count "
          << elem_parts.size() << " != "
          << num_elem_blk ;
    }
    throw std::runtime_error( msg.str() );
  }

  //--------------------------------------------------------------------
  // Read and broadcast element block sizes

  std::vector<int> elem_blk_counts( elem_parts.size() + 1 );

  if ( ! exo_error ) {

    if ( reader ) {

      for ( int i_blk = 0 ; ! exo_error &&
                            i_blk < num_elem_blk ; ++i_blk ) {
        const FilePart & part = * elem_parts[i_blk] ;
        int num_elem_this_blk ;
        int num_node_per_elem ;
        int num_attr_per_elem ;
        char elem_type[ Maximum_Name_Length ];

        exo_func = "ex_get_elem_block" ;
        exo_error = ex_get_elem_block( m_exo_id , part.m_identifier ,
                                       elem_type ,
                                       & num_elem_this_blk ,
                                       & num_node_per_elem ,
                                       & num_attr_per_elem );

        elem_blk_counts[i_blk] = num_elem_this_blk ;
      }
    }

    elem_blk_counts[ num_elem_blk ] = exo_error ;

    broadcast( p_comm, p_read, & elem_blk_counts[0], elem_blk_counts.size() );

    exo_error = elem_blk_counts[ num_elem_blk ];
  }

  //--------------------------------------------------------------------
  // Read and scatter element identifiers and connectivity.
  // Element data is partitioned by element block.
  // Determine size of element data and upper bound size of
  // the needed node map.

  // { elem-index , elem-identifier , { node-index } }

  unsigned size_elem_data_local = 0 ;
  unsigned size_needed_local_nodes = 0 ;
  
  if ( ! exo_error ) { // Sizing pass

    unsigned i = 0 ;

    for ( unsigned i_blk = 0 ; i_blk < elem_parts.size() ; ++i_blk ) {
      const FilePart & part = * elem_parts[i_blk] ;

      const unsigned num_elem_this_blk = elem_blk_counts[i_blk] ;
      const unsigned num_value = 2 + part.m_number_nodes ;
      const unsigned num_items =
        index_processor_count( p_rank, p_size,
                               i, num_elem_this_blk, num_elems_global );
      i += num_elem_this_blk ;

      size_elem_data_local += num_value * num_items ;
      size_needed_local_nodes += part.m_number_nodes * num_items ;
    }
  }

  std::vector<int> elem_data_local( size_elem_data_local );
  std::vector<int> needed_local_nodes( size_needed_local_nodes );

  size_elem_data_local = 0 ;
  size_needed_local_nodes = 0 ;

  if ( ! exo_error ) {

    std::vector<int> data ;
    std::vector<int> connect ;
    std::vector<int> identifiers ;
    std::vector<unsigned> send_size ;

    if ( reader ) { send_size.resize( p_size ); }

    // Load element identifiers and
    // connectivity node indices (not identifiers)

    int i = 0 ;

    for ( unsigned i_blk = 0 ; i_blk < elem_parts.size() ; ++i_blk ) {
      const FilePart & part = * elem_parts[i_blk] ;

      const unsigned num_elem_this_blk = elem_blk_counts[i_blk] ;
      const unsigned num_value = 2 + part.m_number_nodes ;
      const unsigned index_part = 1 + i ;

      unsigned num_item_per_chunk =
        m_max_buffer / ( 2 * num_value * sizeof(int) );

      if ( num_elem_this_blk < num_item_per_chunk ) {
        num_item_per_chunk = num_elem_this_blk ;
      }

      if ( data.size() < num_item_per_chunk * num_value ) {
        data.resize( num_item_per_chunk * num_value );
      }

      if ( reader ) {
        if ( connect.size() < num_item_per_chunk * part.m_number_nodes ) {
          connect.resize( num_item_per_chunk * part.m_number_nodes );
        }
        if ( has_elem_identifiers && identifiers.size() < num_item_per_chunk ) {
          identifiers.resize( num_item_per_chunk );
        }
      }

      const int i_end = i + num_elem_this_blk ;

      while ( i < i_end ) {

        const int index_beg = 1 + i ;

        const int number = num_elem_this_blk < num_item_per_chunk ?
                           num_elem_this_blk : num_item_per_chunk ;

        const unsigned recv_count = num_value *
          index_processor_count( p_rank, p_size, i, number, num_elems_global );

        const unsigned recv_size = sizeof(int) * recv_count ;

        if ( reader ) {

          const unsigned i_last = i + number - 1 ;
          const unsigned p_beg = index_processor_map(p_size,i,num_elems_global);
          const unsigned p_end =
            1 + index_processor_map(p_size,i_last,num_elems_global);

          for ( unsigned p = 0 ;     p < p_beg ;  ++p ) { send_size[p] = 0 ; }
          for ( unsigned p = p_end ; p < p_size ; ++p ) { send_size[p] = 0 ; }

          for ( unsigned p = p_beg ; p < p_end ; ++p ) {
            send_size[p] = sizeof(int) * num_value *
              index_processor_count( p, p_size, i, number, num_elems_global );
          }

          // Read a 'chunk' worth of element identifiers and connectivity
          // Pack into the 'data' array for scattering.

          if ( 0 == exo_error ) {
            const int index_beg_part = 1 + index_beg - index_part ;
            exo_func = "ne_get_n_elem_conn" ;
            exo_error = ne_get_n_elem_conn( m_exo_id , part.m_identifier ,
                                            index_beg_part , number ,
                                            & connect[0] );
          }

          if ( 0 == exo_error && has_elem_identifiers ) {
            exo_func = "ne_get_n_elem_num_map" ; 
            exo_error = ne_get_n_elem_num_map( m_exo_id ,
                                               index_beg , number ,
                                               & identifiers[0] );
          }

          // Copy into scatter data array

          for ( int k = 0 ; k < number ; ++k ) {
            int * const elem_data = & data[ k * num_value ];
            const int * const conn_data = & connect[ k * part.m_number_nodes ];
            elem_data[0] = index_beg + k ;      // Element index
            elem_data[1] = has_elem_identifiers
                           ? identifiers[k] : index_beg + k ; // Identifier
            for ( int j = 0 ; j < part.m_number_nodes ; ++j ) {
              elem_data[2+j] = conn_data[ j ];  // Element-node index
            }
          }
        }

        const unsigned * const p_send_size = reader ? & send_size[0] : NULL ;

        scatter( p_comm , p_read , p_send_size , recv_size , & data[0] );

        for ( unsigned k = 0 ; k < recv_count ; ) {
          elem_data_local[ size_elem_data_local++ ] = data[k++] ; // Index
          elem_data_local[ size_elem_data_local++ ] = data[k++] ; // Id
          for ( int j = 0 ; j < part.m_number_nodes ; ++j ) {
            const int node_index = data[k++] ;
            elem_data_local[ size_elem_data_local++ ] = node_index ;
            needed_local_nodes[ size_needed_local_nodes++ ] = node_index ;
          }
        }

        i += number ;
      }
    }

    { // Sort and unique the index node map, used to request node data.
      std::sort( needed_local_nodes.begin() , needed_local_nodes.end() );
      std::vector<int>::iterator iter =
        std::unique( needed_local_nodes.begin() , needed_local_nodes.end() );
      needed_local_nodes.erase( iter , needed_local_nodes.end() );
    }

    broadcast( p_comm, p_read, & exo_error, 1 );
  }

  //--------------------------------------------------------------------
  // Read and communicate node identifiers and coordinates as needed.

  std::vector<NodeData> node_data_local( needed_local_nodes.size() );

  if ( ! exo_error ) {

    const int num_items_per_chunk =
      m_max_buffer / ( 2 * NodeData::size_of() );

    std::vector<unsigned> send_size ;
    std::vector<int>      ident ;
    std::vector<double>   coord ;
    std::vector<NodeData> data ;

    if ( reader ) {
      send_size.resize( p_size );
      ident.resize( num_items_per_chunk );
      coord.resize( 3 * num_items_per_chunk );
    }

    unsigned count_node_data_local = 0 ;

    std::vector<int>::iterator iter_needed_node_beg =
      needed_local_nodes.begin();

    for ( int i = 0 ; i < num_nodes_global ; ) {
      const int index_beg = 1 + i ;
      int number = num_nodes_global - i ;
      if ( num_items_per_chunk < number ) { number = num_items_per_chunk ; }
      const int index_end = index_beg + number ;

      // Request node data in the span [i,i+number) from the reader.

     iter_needed_node_beg =
        std::lower_bound( iter_needed_node_beg ,
                          needed_local_nodes.end() , index_beg );

      std::vector<int>::iterator iter_needed_node_end =
        std::lower_bound( iter_needed_node_beg ,
                          needed_local_nodes.end() , index_end );

      const unsigned num_needed_node = std::distance( iter_needed_node_beg ,
                                                      iter_needed_node_end );

      const unsigned recv_size = num_needed_node * NodeData::size_of();

      CommGather node_request( p_comm, p_read, sizeof(int) * num_needed_node );

      node_request.send_buffer().pack<int>( & *iter_needed_node_beg ,
                                            num_needed_node );

      iter_needed_node_beg = iter_needed_node_end ;

      node_request.communicate();

      if ( reader ) {

        int    * const m = & ident[ 0 ] ;
        double * const x = & coord[ 0 ] ;
        double * const y = & coord[ number ] ;
        double * const z = & coord[ number * 2 ] ;

        exo_func = "ne_get_n_coord" ;
        exo_error = ne_get_n_coord( m_exo_id, index_beg, number, x, y, z );

        if ( ! exo_error && has_node_identifiers ) {
          exo_func = "ne_get_n_node_num_map" ;
          exo_error = ne_get_n_node_num_map( m_exo_id, index_beg, number, m );
        }

        unsigned total_num_send = 0 ;
        for ( unsigned p = 0 ; p < p_size ; ++p ) {
          CommBuffer & buf_request = node_request.recv_buffer(p);
          const unsigned num_send = buf_request.remaining() / sizeof(int);
          send_size[p] = num_send * NodeData::size_of();
          total_num_send += num_send ;
        }

        if ( data.size() < total_num_send ) { data.resize( total_num_send ); }

        total_num_send = 0 ;

        for ( unsigned p = 0 ; p < p_size ; ++p ) {
          CommBuffer & buf_request = node_request.recv_buffer(p);
          while ( buf_request.remaining() ) {
            int index_node_request ;
            buf_request.unpack<int>( index_node_request );
            const unsigned offset = index_node_request - index_beg ;

            data[ total_num_send ].coord[0] = x[ offset ];
            data[ total_num_send ].coord[1] = y[ offset ];
            data[ total_num_send ].coord[2] = z[ offset ];
            data[ total_num_send ].ident =
              has_node_identifiers ? m[offset] : index_node_request ;
            data[ total_num_send ].index = index_node_request ;
            ++total_num_send ;
          }
        }
      }
      else if ( data.size() < num_needed_node ) {
        data.resize( num_needed_node );
      }

      const unsigned * const p_send_size = reader ? & send_size[0] : NULL ;

      scatter( p_comm , p_read , p_send_size , recv_size , & data[0] );

      for ( unsigned k = 0 ; 0 == exo_error && k < num_needed_node ; ++k ) {
        node_data_local[ count_node_data_local ].coord[0] = data[k].coord[0] ;
        node_data_local[ count_node_data_local ].coord[1] = data[k].coord[1] ;
        node_data_local[ count_node_data_local ].coord[2] = data[k].coord[2] ;
        node_data_local[ count_node_data_local ].ident = data[k].ident ;
        node_data_local[ count_node_data_local ].index = data[k].index ;
        ++count_node_data_local ;
      }

      i += number ;
    }

    broadcast( p_comm , p_read , & exo_error , 1 );
  }

  if ( ! exo_error ) {

    PartSet entity_parts(2);
    entity_parts[0] = & owns_part ;

    // Now have all needed data to create nodes and elements
    // std::vector<int> elem_data_local ;
    // std::vector<NodeData> node_data_local ;

    size_elem_data_local = 0 ;
    unsigned i = 0 ;

    for ( unsigned i_blk = 0 ; i_blk < elem_parts.size() ; ++i_blk ) {
      const FilePart & part = * elem_parts[i_blk] ;

      entity_parts[1] = & part.m_part ;

      const unsigned num_elem_this_blk = elem_blk_counts[i_blk] ;
      const unsigned num_items =
        index_processor_count( p_rank, p_size,
                               i, num_elem_this_blk, num_elems_global );
      i += num_elem_this_blk ;

      for ( unsigned k = 0 ; k < num_items ; ++k ) {

        const int elem_index = elem_data_local[ size_elem_data_local++ ];
        const entity_id_type elem_ident =
            elem_data_local[ size_elem_data_local++ ] ;
        const entity_key_type elem_key = entity_key( Element , elem_ident );

        Entity & elem = M.declare_entity(elem_key,entity_parts,(int)p_rank);

        * elem.data( FS.m_field_elem_index ) = elem_index ;

        for ( int j = 0 ; j < part.m_number_nodes ; ++j ) {
          const int node_index = elem_data_local[ size_elem_data_local++ ];

          // Find the node data

          std::vector<NodeData>::iterator iter_node_data =
            std::lower_bound( node_data_local.begin() ,
                              node_data_local.end() ,
                              node_index ,
                              less_NodeData() );

          if ( iter_node_data == node_data_local.end() ||
               iter_node_data->index != node_index ) {
            std::ostringstream msg ;
            msg << method ;
            msg << " : FAILED to find node_index = " ;
            msg << node_index ;
            throw std::logic_error( msg.str() );
          }

          const entity_id_type node_ident = iter_node_data->ident ;
          const entity_key_type node_key = entity_key( Node , node_ident );

          Entity & node = M.declare_entity(node_key,entity_parts,(int)p_rank);

          M.declare_connection( elem , node , j , method );

          * node.data( FS.m_field_node_index ) = node_index ;

          double * const node_coord = node.data( FS.m_field_node_coord );
          node_coord[0] = iter_node_data->coord[0] ;
          node_coord[1] = iter_node_data->coord[1] ;
          node_coord[2] = iter_node_data->coord[2] ;
        }
      }
    }
  }

  // Nodes and elements are created, discover sharing and generate aura

  comm_mesh_discover_sharing( M );
  comm_mesh_regenerate_aura( M );

  if ( exo_error ) {
    std::ostringstream msg ;
    if ( reader ) {
      msg << method << " FAILED calling " << exo_func << " with " << exo_error ;
    }
    throw std::runtime_error( msg.str() );
  }
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

namespace phdmesh {
namespace exodus {

FileSchema::FileSchema( Schema & arg_schema ,
                        const Field<double,1> & arg_node_coordinates ,
                        const Field<double,1> & arg_elem_attributes ,
                        const std::string     & ,
                        ParallelMachine         ,
                        const unsigned          arg_reader_rank )
  : m_schema( arg_schema ),
    m_io_rank( arg_reader_rank ),
    m_dimension( arg_node_coordinates.max_length() ),
    m_field_node_coord( arg_node_coordinates ),
    m_field_elem_attr(  arg_elem_attributes ),
    m_field_node_index( exo_index( arg_schema , Node ) ),
    m_field_edge_index( exo_index( arg_schema , Edge ) ),
    m_field_face_index( exo_index( arg_schema , Face ) ),
    m_field_elem_index( exo_index( arg_schema , Element ) )
{
}

FileOutput::~FileOutput() {}

FileOutput::FileOutput(
  const FileSchema  & arg_schema ,
  const Mesh        & arg_mesh ,
  const std::string & ,
  const std::string & ,
  const bool ,
  const std::vector<const Field<void,0> * > & ,
  const int * const )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
  {}

void FileOutput::write( double ) {}

FileInput::~FileInput() { }

FileInput::FileInput(
  const FileSchema  & arg_schema ,
        Mesh        & arg_mesh ,
  const std::string & ,
  const std::vector< const Field<void,0> * > & )
  : m_schema( arg_schema ),
    m_mesh( arg_mesh ),
    m_exo_id( 0 ),
    m_counter( 0 ),
    m_max_buffer( 0x0200000 )
{}

}
}

#endif


