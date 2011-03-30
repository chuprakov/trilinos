#ifndef stk_percept_PerceptMesh_hpp
#define stk_percept_PerceptMesh_hpp

#define SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

#include <iostream>
#include <stdexcept>
#include <string>
#include <set>

#include <stk_percept/function/Function.hpp>
#include <stk_percept/Name.hpp>

#include "ShardsInterfaceTable.hpp"

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
//#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <stk_io/IossBridge.hpp>

#include <Intrepid_FieldContainer.hpp>


#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>

#include "Teuchos_RCP.hpp"

#include "PerceptMeshReadWrite.hpp"
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>

#include <stk_percept/SameRankRelation.hpp>


//using namespace Intrepid;

using namespace shards;

namespace Intrepid {
  template<class Scalar, class ArrayScalar>
  class Basis;
}

namespace stk {
  namespace percept {

    typedef mesh::Field<double>                          ScalarFieldType ;
    typedef mesh::Field<double, stk::mesh::Cartesian>    VectorFieldType ;


    using namespace interface_table;

    class PerceptMesh
    {
    public:
      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;
      typedef std::map<unsigned, BasisTypeRCP > BasisTableMap;

      static std::string s_omit_part;

      class GMeshSpec : public Name
      {
      public:
        explicit GMeshSpec(const std::string& name) : Name(name) {}
      };

      struct FieldCreateOrder
      {
        const std::string m_name;
        const unsigned m_entity_rank;
        const std::vector<int> m_dimensions;
        const mesh::Part* m_part;
        FieldCreateOrder();
        FieldCreateOrder(const std::string name, const unsigned entity_rank,
                        const std::vector<int> dimensions, const mesh::Part* part);
      };
      typedef std::vector<PerceptMesh::FieldCreateOrder> FieldCreateOrderVec;


    public:

      //========================================================================================================================
      /// high-level interface

      // ctor constructor
      /// Create a Mesh object that owns its constituent MetaData and BulkData (which are created by this object)
      PerceptMesh(size_t spatialDimension, stk::ParallelMachine comm =  MPI_COMM_WORLD );

      /// reads and commits mesh, editing disabled
      void
      openReadOnly(const std::string& in_filename);

      /// reads but doesn't commit mesh, enabling edit
      void
      open(const std::string& in_filename);

      /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec, Read Only mode, no edits allowed
      void
      newMeshReadOnly(const GMeshSpec gmesh_spec);

      /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec
      void
      newMesh(const GMeshSpec gmesh_spec);

      /// add a field to the mesh
      stk::mesh::FieldBase *
      addField(const std::string& name, const unsigned entity_rank, int vectorDimension=0, const std::string part_name="universal_part");

      stk::mesh::FieldBase *
      getField(const std::string& name);

      /// commits mesh  - any operations done on a non-committed mesh, except to add fields will throw an exception
      void
      commit();

      /// reopens the mesh for editing - warning, this operation writes the mesh to a temp file then re-reads it and
      /// thus recreates the internal MetaData and BulkData
      void
      reopen(const std::string temp_file_name="percept_tmp.e");

      /// commits mesh if not committed and saves it in new file
      void
      saveAs(const std::string& out_filename );

      /// closes this mesh to further changes
      void
      close();

      /// print number of parts and fields, and info on each
      void
      printInfo(std::string header="", int print_level = 0);

      void
      printFields(std::string header="");

      int
      getSpatialDim();

      int
      getNumberElements();

      //========================================================================================================================
      /// low-level interfaces
      /// Create a Mesh object thatPerceptMesh* mesh_data doesn't own its constituent MetaData and BulkData, pointers to which are adopted
      /// by this constructor.
      PerceptMesh(const stk::mesh::MetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted=true);

      ~PerceptMesh() ;
      void init ( stk::ParallelMachine comm);
      void destroy();

      /// reads the given file into a temporary model and prints info about it
      void dump(const std::string& file="");
      void dumpElements(const std::string& partName = "");

      unsigned getRank() { return getBulkData()->parallel_rank(); }
      unsigned getParallelSize() { return getBulkData()->parallel_size(); }
      bool isGhostElement(const stk::mesh::Entity& element)
      {
        //throw std::runtime_error("not impl"); // FIXME
        bool isGhost = element.owner_rank() != getRank();
        return isGhost;
      }

      stk::mesh::Entity & createOrGetNode(stk::mesh::EntityId nid, double* x=0);

      void createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity *>& requested_entities);

      const mesh::Part*
      getPart(const std::string& part_name) ;

      mesh::Part*
      getNonConstPart(const std::string& part_name);

      static double * field_data(const stk::mesh::FieldBase *field, const stk::mesh::Bucket & bucket, unsigned *stride=0);
      static double * field_data(const stk::mesh::FieldBase *field, const mesh::Entity& node, unsigned *stride=0);

      static inline double *
      field_data_inlined(const FieldBase *field, const mesh::Entity& node)
      {
        return field_data(field, node);
//         return 
//           field->rank() == 0 ? 
//           stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , node )
//           :
//           stk::mesh::field_data( *static_cast<const VectorFieldType *>(field) , node );
      }


      double * node_field_data(stk::mesh::FieldBase *field, const mesh::EntityId node_id);

      stk::mesh::BulkData * getBulkData();
      stk::mesh::MetaData * getMetaData();

      static BasisTypeRCP getBasis(shards::CellTopology& topo);
      static void setupBasisTable();

      void nodalOpLoop(GenericFunction& nodalOp, stk::mesh::FieldBase *field);
      void elementOpLoop(ElementOp& elementOp, stk::mesh::FieldBase *field=0, stk::mesh::Part *part = 0);
      void bucketOpLoop(BucketOp& bucketOp, stk::mesh::FieldBase *field=0, stk::mesh::Part *part = 0);

      static unsigned size1(const stk::mesh::Bucket& bucket) { return bucket.size(); }
      static unsigned size1(const stk::mesh::Entity& element) { return 1; }

      /// \brief Fill the array cellNodes(numCells, numNodesPerCell, nDof) with DOF values from the given Field
      /// The stride of the data (last dimension in cellNodes) is taken to be that of the field's stride; however,
      /// in some cases you may want to pass in an array where nDof is less than the stride (e.g., pull out 2
      /// coordinates from a 3D coordinate field).  In that case, the dataStride argument can be set (to e.g. "2").
      template<class ArrayType>
      static void fillCellNodes( const stk::mesh::Bucket &bucket,
                                 //stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                 //VectorFieldType& coord_field,
                                 FieldBase* field,
                                 ArrayType& cellNodes, unsigned dataStride=0 );

      /// \brief see comment for fillCellNodes(Bucket& ...)
      template<class ArrayType>
      static void fillCellNodes( const stk::mesh::Entity &element,
                                 //stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                 //VectorFieldType& coord_field,
                                 FieldBase* field,
                                 ArrayType& cellNodes, unsigned dataStride=0 );

      static void findMinMaxEdgeLength(const mesh::Bucket &bucket,  stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                       FieldContainer<double>& elem_min_edge_length, FieldContainer<double>& elem_max_edge_length);

      VectorFieldType* getCoordinatesField() {
        // this should have been set by a previous internal call to setCoordinatesField
        return m_coordinatesField;
      }

      static void
      element_side_nodes( const Entity & elem , int local_side_id, EntityRank side_entity_rank, std::vector<Entity *>& side_node_entities );

      static void
      element_side_permutation(const Entity& element, const Entity& side, unsigned iSubDimOrd, int& returnedIndex, int& returnedPolarity);

      // FIXME
      SameRankRelation& adapt_parent_to_child_relations() { return m_adapt_parent_to_child_relations; }

      bool
      isBoundarySurface(mesh::Part& block, mesh::Part& surface);

      /// here @param thing is a Part, Bucket, Entity, or Field or BulkData
      template<class T>
      static
      const stk::mesh::fem::FEMMetaData& get_fem_meta_data(const T& thing) 
      { 
        const stk::mesh::MetaData& meta = MetaData::get(thing);
        const stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get ( meta );
        return fem_meta;
      }

      /// here @param thing is a Part, Bucket, Entity
      template<class T>
      static
      const CellTopologyData * const get_cell_topology(const T& thing) 
      { 
        const CellTopologyData * const cell_topo_data = fem::get_cell_topology_new(thing).getCellTopologyData();
        return cell_topo_data;
      }



    private:

      /// reads meta data, commits it, reads bulk data
      void readModel( const std::string& in_filename );

      /// read with no commit
      void readMetaDataNoCommit( const std::string& in_filename );

      /// create with no commit
      void createMetaDataNoCommit( const std::string& gmesh_spec);

      void commitMetaData();

      /// read the bulk data (no op in create mode)
      void readBulkData();

      /// Convenience method to read a model's meta data, create some new fields, commit meta data then read the bulk data
      /// deprecated
      void readModelAndCreateOptionalFields(const std::string file, bool print,  FieldCreateOrderVec create_field);

      /// after the meta data is read or created, create some fields using this method, then you can commit and read bulk data(if in
      /// reading mode, else no need to read bulk data in create mode)
      // deprecated
      void createFields(bool print, FieldCreateOrderVec create_field = FieldCreateOrderVec());

      /// Cache internal pointer to coordinate field
      void setCoordinatesField();

      // write in exodus format to given file
      void writeModel( const std::string& out_filename );

      stk::mesh::FieldBase * createField(const std::string& name, const unsigned entity_rank, const std::vector<int>& dimensions,
                                         const stk::mesh::Part* arg_part=0);

      //static void transformMesh(GenericFunction& coordinate_transform);

    private:
      stk::mesh::fem::FEMMetaData *         m_femMetaData;
      stk::mesh::MetaData *                 m_metaData;
      stk::mesh::BulkData *                 m_bulkData;
      stk::io::util::Gmesh_STKmesh_Fixture* m_fixture;
      Ioss::Region *                        m_iossRegion;
      VectorFieldType*                      m_coordinatesField;
      int                                   m_spatialDim;
      bool                                  m_ownData;
      bool                                  m_isCommitted;
      bool                                  m_isOpen;
      bool                                  m_isInitialized;
      bool                                  m_isAdopted;
      bool                                  m_dontCheckState;
      std::string                           m_filename;
      stk::ParallelMachine                  m_comm;

      //static std::map<unsigned, BasisType *> m_basisTable;
      static BasisTableMap m_basisTable;

      SameRankRelation m_adapt_parent_to_child_relations;

      void checkStateSpec(const std::string& function, bool cond1=true, bool cond2=true, bool cond3=true);

      void checkState(const std::string& function) {
        return checkStateSpec(function, m_isOpen, m_isInitialized, m_isCommitted);
      }

    }; // class PerceptMesh


    template<>
    const CellTopologyData * const 
    PerceptMesh::get_cell_topology(const Part& part) ;



#if 0
    inline
    std::string &operator<<(std::string& out, const char *str)
    {
      return out.append(str);
    }
#endif

    // static
    template<class ArrayType>
    void PerceptMesh::fillCellNodes( const mesh::Bucket &bucket,
                                  //stk::mesh::Cartesian>& coord_field,
                                  //VectorFieldType& coord_field,
                                  FieldBase* field,
                                  ArrayType& cellNodes, unsigned dataStrideArg)
    {
      unsigned number_elems = bucket.size();
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
      const CellTopologyData * const bucket_cell_topo_data = PerceptMesh::get_cell_topology(bucket);
#else
      const CellTopologyData * const bucket_cell_topo_data = stk::mesh::fem::get_cell_topology(bucket).getCellTopologyData();
#endif

      CellTopology cell_topo(bucket_cell_topo_data);
      //unsigned numCells = number_elems;
      unsigned numNodes = cell_topo.getNodeCount();
      //unsigned spaceDim = cell_topo.getDimension();

      unsigned dataStride = dataStrideArg;
      if (!dataStrideArg)
        {
          const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::Node, MetaData::get(*field).universal_part());
          dataStride = r.stride[0] ;
        }
      //std::cout << "bucket dataStride= " << dataStride << std::endl;

      for ( unsigned iElemInBucketOrd = 0 ; iElemInBucketOrd < number_elems ; ++iElemInBucketOrd)
        {
          mesh::Entity & elem = bucket[iElemInBucketOrd] ;

          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          const mesh::PairIterRelation elem_nodes = elem.relations( mesh::Node );

          // FIXME: fill field data (node coordinates)
          for (unsigned iNodeOrd = 0; iNodeOrd < numNodes; iNodeOrd++)
            {
              mesh::Entity& node = *elem_nodes[iNodeOrd].entity();
              double * node_coord_data = PerceptMesh::field_data( field , node);

              for (unsigned iDOFOrd = 0; iDOFOrd < dataStride; iDOFOrd++)
                {
                  cellNodes(iElemInBucketOrd, iNodeOrd, iDOFOrd) = node_coord_data[iDOFOrd];
                }
            }


        }

    }

    // static
    template<class ArrayType>
    void PerceptMesh::fillCellNodes( const stk::mesh::Entity &element,
                                  //stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                  //VectorFieldType& coord_field,
                                  FieldBase* field,
                                  ArrayType& cellNodes,
                                  unsigned dataStrideArg )
    {
      unsigned dataStride = dataStrideArg;
      if (!dataStrideArg)
        {
          const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::Node, MetaData::get(*field).universal_part());
          dataStride = r.stride[0] ;
        }
      //std::cout << "element dataStride= " << dataStride << std::endl;
      const mesh::PairIterRelation element_nodes = element.relations( mesh::Node );
      unsigned numNodes = element_nodes.size();

      unsigned iCell = 0;
      for (unsigned iNodeOrd = 0; iNodeOrd < numNodes; iNodeOrd++)
        {
          mesh::Entity& node = *element_nodes[iNodeOrd].entity();
          double * node_coord_data = PerceptMesh::field_data( field , node);

          for (unsigned iDOFOrd = 0; iDOFOrd < dataStride; iDOFOrd++)
            {
              cellNodes(iCell, iNodeOrd, iDOFOrd) = node_coord_data[iDOFOrd];
            }
        }
    }


  }
}

#endif
