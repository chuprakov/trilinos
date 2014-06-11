/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_FEINTERPOLATE_HPP
#define  STK_FEINTERPOLATE_HPP

#include <boost/shared_ptr.hpp>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/environment/ReportHandler.hpp>


namespace stk {
namespace transfer {

template <class MESHA, class MESHB> class FEInterpolate {

public :

  typedef MESHA                            MeshA;
  typedef MESHB                            MeshB;
  typedef typename MeshA::EntityKey        EntityKeyA;
  typedef typename MeshB::EntityKey        EntityKeyB;
  typedef typename MeshA::EntityProc       EntityProcA;
  typedef typename MeshB::EntityProc       EntityProcB;

  typedef std::pair<EntityProcB, EntityProcA>             EntityProcRelation;
  typedef std::vector<EntityProcRelation>                 EntityProcRelationVec;

  typedef std::multimap<EntityKeyB, EntityKeyA>        EntityKeyMap;

  enum { Dimension = MeshA::Dimension };
  
  static void filter_to_nearest(EntityKeyMap    &BtoA,
                                const MeshA     &FromElem,
                                      MeshB     &ToPoints);
  
  static void apply (MeshB               &meshb,
                     const MeshA         &mesha,
                     const EntityKeyMap &RangeToDomain);
  

  struct Record : public MeshB::Record {
    virtual ~Record(){}
    std::vector<double> P;
  };

protected :

  static typename EntityKeyMap::iterator determine_best_fit(
                                         const typename EntityKeyMap::iterator begin,
                                         const typename EntityKeyMap::iterator end,
                                         const EntityKeyB          current,
                                               MeshB             &ToPoints,
                                         const MeshA             &FromElem);
  
private :

  enum { dim_eq = StaticAssert<static_cast<unsigned>(MeshB::Dimension)==static_cast<unsigned>(MeshA::Dimension)>::OK };
  enum { dim_3  = StaticAssert<               3==MeshA::Dimension>::OK };
  
};

template <class MESHA, class MESHB> typename FEInterpolate<MESHA,MESHB>::EntityKeyMap::iterator 
FEInterpolate<MESHA,MESHB>::determine_best_fit(
                               const typename EntityKeyMap::iterator begin,
                               const typename EntityKeyMap::iterator end,
                               const EntityKeyB          current,
                                     MeshB             &ToPoints,
                               const MeshA             &FromElems) {


  typename EntityKeyMap::iterator nearest = end;
  double min_distance     = std::numeric_limits<double>::max();
  const double *to_coords = ToPoints.coord(current);

  for (typename EntityKeyMap::iterator d=begin; d != end; ++d) {
    std::vector<double> parametric;  
    const EntityKeyA domain_index = d->second;
    const double dist = FromElems.parametric_coord(parametric, to_coords, domain_index);
    if (dist < min_distance) {
      min_distance = dist;
      Record *record = ToPoints.template database<Record>(d->first); 
      record->P=parametric;
      nearest = d;
    }
  }
  return nearest;
}

template <class MESHA, class MESHB> void FEInterpolate<MESHA,MESHB>::filter_to_nearest(
                                    EntityKeyMap  &BtoA,
                                    const MeshA   &mesha,
                                          MeshB   &meshb) {
  typedef typename EntityKeyMap::iterator iterator;
  for (iterator j=BtoA.begin(); j!=BtoA.end(); ) {
    std::pair<iterator, iterator> keys=BtoA.equal_range(j->first);
    typename EntityKeyMap::iterator n = determine_best_fit(keys.first, keys.second, j->first, meshb, mesha);
    if (n != keys.first ) BtoA.erase(keys.first, n);
    if (n != keys.second) BtoA.erase(++n, keys.second);
    j = keys.second;
  }
}


template <class MESHA, class MESHB>  void FEInterpolate<MESHA,MESHB>::apply (
                          MeshB        &meshb,
                    const MeshA        &mesha,
                    const EntityKeyMap &RangeToDomain){

  typedef typename EntityKeyMap::const_iterator const_iterator;

  const unsigned numValsa = mesha.num_values();
  const unsigned numValsb = meshb.num_values();
    ThrowRequireMsg (numValsb == numValsa,  
      __FILE__<<":"<<__LINE__<<" Found "<<numValsa<<" values for mesh a and "<<numValsb<<" for mesh b."
      <<" These should be the same.");

  for (const_iterator j,i=RangeToDomain.begin(); i!=RangeToDomain.end(); i=j)  {
    j = i;
    while (j != RangeToDomain.end() && i->first == j->first) ++j;
    const unsigned num_relations = distance(i, j);

    const EntityKeyB to_key = i->first;
    ThrowRequireMsg (1 == num_relations,  
      __FILE__<<":"<<__LINE__<<" Expected "<<1<<" relation."<<" Found:"<<num_relations<<" for Key:"<<to_key);
    const EntityKeyA domain_index         = i->second;


    const Record *record = meshb.template database<Record>(to_key); 
    const std::vector<double> &parametric = record->P;

    std::vector<std::vector<double> > val;
    mesha.eval_parametric(val, parametric, domain_index);

    for (unsigned f=0; f<numValsb; ++f)  {
      const unsigned   to_field_size =   meshb.value_size(to_key,       f);
      const unsigned from_field_size =   mesha.value_size(domain_index, f);
      const unsigned field_size = std::min(to_field_size, from_field_size);

      double  *c = meshb.value(to_key, f);
      for (unsigned n=0; n<field_size; ++n) c[n] = val[f][n];
    }
  }
}

} }
#endif

