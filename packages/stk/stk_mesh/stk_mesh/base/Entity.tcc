/** \addtogroup stk_mesh_module
 * \{
 */

namespace stk {
namespace adapt {
struct my_tuple_hash;
}
}

#ifdef SIERRA_MIGRATION

namespace stk {
namespace mesh {
union Entity;
class BulkData;
}
}

namespace sierra {
namespace Fmwk {

class MeshObjRoster;
class MeshObjSharedAttr;
class MeshBulkData;

extern const unsigned int INVALID_LOCAL_ID;
extern const stk::mesh::RelationIterator INVALID_RELATION_ITR;

namespace detail {
bool set_attributes( MeshBulkData& meshbulk, stk::mesh::Entity , const int , const MeshObjSharedAttr*, const int);
bool set_attributes( MeshBulkData& meshbulk, stk::mesh::Entity , const MeshObjSharedAttr*, const int);
bool update_relation( stk::mesh::Entity, const stk::mesh::RelationIterator ir, const bool back_rel_flag, MeshBulkData& bulk);
}

namespace roster_only {
void destroy_meshobj(stk::mesh::Entity, MeshBulkData& meshbulk );
}

const MeshObjSharedAttr * get_shared_attr(const stk::mesh::Entity mesh_obj, const stk::mesh::BulkData& meshbulk);
bool insert_relation( stk::mesh::Entity , const stk::mesh::RelationType,  stk::mesh::Entity , const unsigned, const unsigned, const bool, MeshBulkData &);
bool remove_relation(stk::mesh::Entity , const stk::mesh::RelationIterator, MeshBulkData &);
}
}
#endif

namespace stk {
namespace mesh {

namespace impl {
class EntityRepository;
}

//----------------------------------------------------------------------
// NKC OPT, is this now just the local_offset?  IS the shift ever still needed
union Entity
{
  enum Entity_t {
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
    BULK_DATA_ID_SHIFT = 48ull,
#endif
    InvalidEntity = 0ull,
    MinEntity = 1ull,
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
    MaxEntity = (1ull << BULK_DATA_ID_SHIFT) - 1,
    ForceEnumToBeUint64 = ~0ull,
#else
    MaxEntity = ~0ull,
#endif
    LOCAL_OFFSET_MASK = MaxEntity
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
    , BULK_DATA_ID_MASK = (1ull << 16) - 1
#endif
  };

  uint64_t m_value;

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  struct {
    unsigned char bulk_data_id_part;
    uint32_t local_offset32_part;
  } view1;

  struct {
    uint32_t local_offset32_part;
    unsigned char bulk_data_id_part;
  } view2;
#endif

  Entity operator=(Entity_t val) { m_value = val; return *this;}

  /** \brief local_offset is this entity's offset into all local entities of the same rank.
   * An entity's local_offset will generally remain unchanged through mesh-modification cycles,
   * which means that the set of local_offsets may not be compact or contiguous if some
   * entities have been deleted. (local_offsets for deleted entities are no longer valid.)
   * Thus, local_offset is not suitable for use as an equation index for linear-system operations.
   * See local_id() below.
   */
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  size_t local_offset() const { return static_cast<size_t>(m_value & LOCAL_OFFSET_MASK); }
#else
  size_t local_offset() const { return static_cast<size_t>(m_value); }
#endif

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  int bulk_data_id() const { return BULK_DATA_ID_MASK & (m_value >> BULK_DATA_ID_SHIFT); }
#endif

  bool is_local_offset_valid() const { return local_offset() > 0; }

  /** This setter should only be called by the BulkData class when creating entities.
   * Erroneous calls will lead to undefined (and probably disastrous) behavior.
   */
  void set_local_offset(size_t localOffset) {
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
    int bdi = bulk_data_id();
#endif
    m_value = static_cast<Entity_t>(localOffset);
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
    set_bulk_data_id(bdi);
#endif
  }

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  void set_bulk_data_id(int bdi) {
    int64_t val = bdi;
    val = val << BULK_DATA_ID_SHIFT;
    m_value |= val;
  }
#endif

  bool operator==(Entity entity) const { return m_value == entity.m_value; }

  bool operator==(Entity_t val) const { return m_value == val; }

  bool operator!=(Entity entity) const { return m_value != entity.m_value; }

  bool operator<(Entity entity) const { return m_value < entity.m_value; }

#ifdef SIERRA_MIGRATION

  enum ObjectTypeEnum {
    NODE       = 0,
    EDGE       = 1,
    FACE       = 2,
    ELEMENT    = 3,
    CONSTRAINT = 4,
    NUM_TYPES  = 5,
    BASE_CLASS = 0x00ff
  };

  static std::string TypeToString (ObjectTypeEnum type);

#endif

  //
  // NEED TO REFACTOR CALLERS TO ELIMINATE THE FOLLOWING
  //

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  bool is_valid() const;
  EntityState state() const;
  EntityRank entity_rank() const;
  EntityId identifier() const;
  const EntityKey key() const;
  Bucket & bucket() const;
  Bucket * bucket_ptr() const;
  size_t bucket_ordinal() const;
  size_t synchronized_count() const;
  int owner_rank() const;
#endif


#ifdef SIERRA_MIGRATION
  friend class sierra::Fmwk::MeshObjRoster;
  // These are free functions to facilitate the stk migration:
  friend bool sierra::Fmwk::detail::set_attributes( sierra::Fmwk::MeshBulkData&, Entity , const int, const sierra::Fmwk::MeshObjSharedAttr *, const int);
  friend bool sierra::Fmwk::detail::set_attributes( sierra::Fmwk::MeshBulkData&, Entity , const sierra::Fmwk::MeshObjSharedAttr *, const int);
  friend bool sierra::Fmwk::insert_relation( Entity , const stk::mesh::RelationType, Entity , const unsigned, const unsigned, const bool, sierra::Fmwk::MeshBulkData &);
  friend bool sierra::Fmwk::detail::update_relation( Entity, const stk::mesh::RelationIterator ir, const bool back_rel_flag, sierra::Fmwk::MeshBulkData& bulk);
  friend bool sierra::Fmwk::remove_relation(Entity , const stk::mesh::RelationIterator, sierra::Fmwk::MeshBulkData &);
  friend void sierra::Fmwk::roster_only::destroy_meshobj(stk::mesh::Entity, sierra::Fmwk::MeshBulkData& meshbulk );
#endif
};

std::ostream & operator << ( std::ostream & , const Entity & );

//
// THIS CLASS STAYS BUT NEEDS TO BE DECLARED AFTER Entity.
//
class EntityEqual
{
public:
  bool operator()(const Entity lhs, const Entity rhs) const {
    return lhs == rhs;
  }
};

/** \} */

} // namespace mesh
} // namespace stk
