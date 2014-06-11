/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_GEOMETRICTRANSFER_HPP
#define  STK_GEOMETRICTRANSFER_HPP

#include <set>
#include <vector>
#include <string>
#include <algorithm>

#include <boost/shared_ptr.hpp>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_transfer/TransferBase.hpp>

namespace stk {
namespace transfer {

template <class T>
typename T::EntityKeyMap::iterator  insert (typename T::EntityKeyMap &map,
                                      const typename T::EntityKeyMap::key_type &k,
                                      const typename T::EntityKeyA &a) {
  const typename T::EntityKeyMap::mapped_type m(a);
  const typename T::EntityKeyMap::value_type  v(k, m);
  const typename T::EntityKeyMap::iterator it = map.insert(v);
  return it;
}

template <class INTERPOLATE> class GeometricTransfer : public TransferBase {

public :
  
  typedef INTERPOLATE                                     InterpolateClass;
  typedef typename InterpolateClass::MeshA                MeshA;
  typedef typename InterpolateClass::MeshB                MeshB;
  typedef typename InterpolateClass::EntityKeyA           EntityKeyA;
  typedef typename InterpolateClass::EntityKeyB           EntityKeyB;
  typedef typename InterpolateClass::EntityKeyMap         EntityKeyMap;

  typedef typename InterpolateClass::EntityProcA          EntityProcA;
  typedef typename InterpolateClass::EntityProcB          EntityProcB;

  typedef typename InterpolateClass::EntityProcRelation       EntityProcRelation;
  typedef typename InterpolateClass::EntityProcRelationVec    EntityProcRelationVec;

  typedef typename std::set     <EntityKeyA>              EntityKeySetA;
  typedef typename std::set     <EntityKeyB>              EntityKeySetB;
  typedef typename MeshA::BoundingBox                     BoundingBoxA;
  typedef typename MeshB::BoundingBox                     BoundingBoxB;
  
  
  enum {Dimension = MeshA::Dimension};
  
  GeometricTransfer(boost::shared_ptr<MeshA> &mesha, 
                    boost::shared_ptr<MeshB> &meshb,
                    const std::string &name,
                    const double expansion_factor = 1.5);
  virtual ~GeometricTransfer(){};
  virtual void coarse_search();
  virtual void communication();
  virtual void local_search();
  virtual void apply();

  void determine_entities_to_copy(typename MeshA::EntityProcVec &entities_to_copy) const;
  const boost::shared_ptr<MeshA> mesha() const {return m_mesha;}
  const boost::shared_ptr<MeshB> meshb() const {return m_meshb;}

protected :
  boost::shared_ptr<MeshA>               m_mesha;
  boost::shared_ptr<MeshB>               m_meshb;

  void copy_domain_to_range_processors();
  void localize_entity_key_map();
  
private :
  
  const std::string     m_name;
  const double          m_expansion_factor;
  
  EntityProcRelationVec m_global_range_to_domain;
  EntityKeyMap          m_local_range_to_domain;
  
  template <bool B, typename T = void> struct Enable_If          { typedef T type; };
  template <        typename T       > struct Enable_If<false, T>{                 };

  template <typename T> struct optional_functions {
  private : 
    template <typename X, X> class check {};
    template <typename Type, Type Ptr> struct MemberHelperClass;

    template <typename X> static long copy_ent(...);
    template <typename X> static char copy_ent(
      check<
        void (X::*)(const typename MeshA::EntityProcVec&, const std::string&), 
        &X::copy_entities
      >*);

    template <typename X> static long pcsf(...);
    template <typename X> static char pcsf(
      check<
        void(*)(EntityProcRelationVec&, const MeshA&, const MeshB&), 
        &X::post_coarse_search_filter
      >*);
  public :
    static const bool copy_entities             = (sizeof(copy_ent<T>(0)) == sizeof(char));
    static const bool post_coarse_search_filter = (sizeof(pcsf    <T>(0)) == sizeof(char));
  };
  
  template <typename T> typename Enable_If<optional_functions<T>::copy_entities>::type
       copy_entities(T                                   &mesh,
                     const typename MeshA::EntityProcVec &entities_to_copy,
                     const std::string                   &transfer_name) const 
  { mesh.copy_entities(entities_to_copy, transfer_name); }

  template <typename T> typename Enable_If<!optional_functions<T>::copy_entities>::type
       copy_entities(T                                   &mesh,
                     const typename MeshA::EntityProcVec &entities_to_copy,
                     const std::string                   &transfer_name) const {
    ThrowErrorMsg(__FILE__<<":"<<__LINE__<<" Error: copy_entities undefinded in this class.");
  }

  template <typename T> static typename Enable_If<optional_functions<T>::post_coarse_search_filter>::type
       post_coarse_search_filter(EntityProcRelationVec &BtoA,
                     const MeshA                       &mesha,
                     const MeshB                       &meshb) {
    T::post_coarse_search_filter(BtoA, mesha, meshb);
  }
  template <typename T> static typename Enable_If<!optional_functions<T>::post_coarse_search_filter>::type
       post_coarse_search_filter(EntityProcRelationVec &BtoA,
                     const MeshA                       &msha,
                     const MeshB                       &meshb) {
    ThrowErrorMsg(__FILE__<<":"<<__LINE__<<" Error: post_coarse_search_filter undefinded in this class.");
  }




  void coarse_search(EntityProcRelationVec   &RangeToDomain,
                     const MeshA             &mesha,
                     const MeshB             &meshb,
                     const double             expansion_factor) const ;

  struct compare {
    bool operator()(const BoundingBoxB &a, const EntityProcB  &b) const;
    bool operator()(const EntityProcB  &a, const BoundingBoxB &b) const;
  };

  void delete_range_points_found(std::vector<BoundingBoxB>            &range_vector,
                                 const EntityProcRelationVec          &del) const ;
  
  enum { dim_eq = StaticAssert<static_cast<unsigned>(MeshB::Dimension)==static_cast<unsigned>(MeshA::Dimension)>::OK };
  
};



template <class INTERPOLATE> GeometricTransfer<INTERPOLATE>::GeometricTransfer (boost::shared_ptr<MeshA> &mesha,
                                                                                boost::shared_ptr<MeshB> &meshb,
                                                                                const std::string        &name,
                                                                                const double              expansion_factor) :
  m_mesha(mesha),
  m_meshb(meshb),
  m_name (name) ,
  m_expansion_factor(expansion_factor) {}


template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::coarse_search() {

  m_global_range_to_domain.clear();
  coarse_search(m_global_range_to_domain, 
                *m_mesha, 
                *m_meshb, 
                m_expansion_factor);
}
template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::communication() {

  ParallelMachine comm = m_mesha->comm();
  const unsigned p_size = parallel_machine_size(comm);
  if (1 < p_size) {
    if (optional_functions<MeshA>::copy_entities) copy_domain_to_range_processors();
    else ThrowRequireMsg (optional_functions<MeshA>::copy_entities,
             __FILE__<<":"<<__LINE__<<" Still working on communicaiton capabilities.");
  }
}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::local_search() {
  localize_entity_key_map();
  INTERPOLATE::filter_to_nearest(m_local_range_to_domain, *m_mesha, *m_meshb); 
}


template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::apply(){
  m_mesha->update_values();
  INTERPOLATE::apply(*m_meshb, *m_mesha, m_local_range_to_domain);
  m_meshb->update_values();
}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::determine_entities_to_copy(
                       typename MeshA::EntityProcVec   &entities_to_copy) const {
  entities_to_copy.clear();

  ParallelMachine comm = m_mesha->comm();
  const unsigned my_rank = parallel_machine_rank(comm);
  
  const typename EntityProcRelationVec::const_iterator end=m_global_range_to_domain.end();
  for (typename EntityProcRelationVec::const_iterator i=m_global_range_to_domain.begin(); i!=end; ++i) {
    const unsigned            domain_owning_rank = i->second.proc;
    const unsigned             range_owning_rank = i->first.proc;
    if (domain_owning_rank == my_rank && range_owning_rank != my_rank) {
      const EntityKeyA entity = i->second.ident;
      const typename MeshA::EntityProc ep(entity, range_owning_rank);
      entities_to_copy.push_back(ep);
    }   
  }
  std::sort(entities_to_copy.begin(), entities_to_copy.end());
  typename MeshA::EntityProcVec::iterator del = std::unique(entities_to_copy.begin(), entities_to_copy.end());
  entities_to_copy.resize(std::distance(entities_to_copy.begin(), del));
}
  

template <class INTERPOLATE> void
GeometricTransfer<INTERPOLATE>::localize_entity_key_map()  {

  ParallelMachine comm = m_mesha->comm();
  const unsigned my_rank = parallel_machine_rank(comm);
  
  m_local_range_to_domain.clear();
  for (typename EntityProcRelationVec::const_iterator i=m_global_range_to_domain.begin(); i!=m_global_range_to_domain.end(); ++i) {
    const unsigned range_owning_rank = i->first.proc;
    if (range_owning_rank == my_rank) insert<INTERPOLATE>(m_local_range_to_domain, i->first.ident, i->second.ident);
  }   
}

template <class INTERPOLATE> void
GeometricTransfer<INTERPOLATE>::copy_domain_to_range_processors()  {
  
  typename MeshA::EntityProcVec entities_to_copy ;

  determine_entities_to_copy(entities_to_copy);
  copy_entities(*m_mesha, entities_to_copy, m_name);
}

template <class INTERPOLATE> bool GeometricTransfer<INTERPOLATE>::compare::operator()(const BoundingBoxB &a, const EntityProcB &b) const {
  return a.key < b;
}
template <class INTERPOLATE> bool GeometricTransfer<INTERPOLATE>::compare::operator()(const EntityProcB &a, const BoundingBoxB &b) const {
  return a < b.key;
}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::delete_range_points_found(
                               std::vector<BoundingBoxB>            &range_vector,
                               const EntityProcRelationVec          &del) const {

  std::vector<EntityProcB> range_entities_found;
  range_entities_found.reserve(del.size());
  for (typename EntityProcRelationVec::const_iterator i=del.begin(); i!=del.end(); ++i) {
    range_entities_found.push_back(i->first);
  }
  {
    std::sort(range_entities_found.begin(), range_entities_found.end());
    const typename std::vector<EntityProcB>::iterator it = std::unique(range_entities_found.begin(), range_entities_found.end());
    range_entities_found.resize(it-range_entities_found.begin());
  }
  
  std::vector<BoundingBoxB> difference(range_vector.size());
  {
    const typename std::vector<BoundingBoxB>::iterator it = 
      std::set_difference(
        range_vector.        begin(), range_vector.        end(),
        range_entities_found.begin(), range_entities_found.end(),
        difference.begin(), compare());
    difference.resize(it-difference.begin());
  }
  swap(difference, range_vector);  
}

template <class INTERPOLATE>  void GeometricTransfer<INTERPOLATE>::coarse_search
(EntityProcRelationVec   &range_to_domain,
 const MeshA             &mesha,
 const MeshB             &meshb,
 const double            expansion_factor) const {
  
  std::vector<BoundingBoxB> range_vector;
  std::vector<BoundingBoxA> domain_vector;

  mesha.bounding_boxes(domain_vector);
  meshb.bounding_boxes(range_vector);

  search::FactoryOrder order;
  order.m_communicator = mesha.comm();

  unsigned range_vector_not_empty = !range_vector.empty();
  stk::all_reduce( mesha.comm(), stk::ReduceSum<1>(&range_vector_not_empty));
  while (range_vector_not_empty) { // Keep going until all range points are processed.
    // Slightly confusing: coarse_search documentation has domain->range
    // relations sorted by domain key.  We want range->domain type relations
    // sorted on range key. It might appear we have the arguments revered
    // in coarse_search call, but really, this is what we want.  
    EntityProcRelationVec rng_to_dom;
    search::coarse_search(rng_to_dom, domain_vector, range_vector, order);

    if (optional_functions<InterpolateClass>::post_coarse_search_filter) 
      post_coarse_search_filter<InterpolateClass>(rng_to_dom, mesha, meshb); 
    range_to_domain.insert(range_to_domain.end(), rng_to_dom.begin(), rng_to_dom.end());

    delete_range_points_found(range_vector, rng_to_dom); 

    for (typename std::vector<BoundingBoxB>::iterator i=range_vector.begin(); i!=range_vector.end(); ++i) {
      // If points were missed, increase search radius.
      i->scale(expansion_factor);
    }
    range_vector_not_empty = !range_vector.empty();
    stk::all_reduce( mesha.comm(), stk::ReduceSum<1>(&range_vector_not_empty));
  } 
  sort (range_to_domain.begin(), range_to_domain.end());
}

}
}


#endif

