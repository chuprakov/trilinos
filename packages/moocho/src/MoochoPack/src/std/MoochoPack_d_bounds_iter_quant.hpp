// ////////////////////////////////////////////////////////////////////
// d_bounds_iter_quant.h

#ifndef D_BOUNDS_ITE_QUANT_HH
#define D_BOUNDS_ITE_QUANT_HH

#include "ActSetStats.h"
#include "GeneralIterationPack/include/CastIQMember.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"

namespace ReducedSpaceSQPPack {

/// Name given to the bounds for d iteration quanity
extern const std::string d_bounds_name;

///
/** Encapsulation class for spare bounds
  */
struct SparseBounds {
	SpVector	l;	// Lower bounds
	SpVector	u;	// Upper bounds
};

///
/** Class for object that attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name act_set_stats_name.
  */
class d_bounds_iq_member
	: public CastIQMember<SparseBounds>
{
public:
    d_bounds_iq_member()
    	: CastIQMember<SparseBounds>(d_bounds_name)
    {}
};

}	// end namespace ReducedSpaceSQPPack

#endif	// D_BOUNDS_ITE_QUANT_HH
