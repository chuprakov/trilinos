// ////////////////////////////////////////////////////////////
// SortByDescendingAbsValue.h

#ifndef SORT_BY_DESCENDING_ABS_VALUE_H
#define SORT_BY_DESCENDING_ABS_VALUE_H

#include <math.h>
#include "SpVectorClass.h"

namespace SparseLinAlgPack {

///
/** Function object class for sorting a sparse vectors in descending order
  * by abs(v(i)).
  */
class SortByDescendingAbsValue {
public:
	bool operator()( const SparseLinAlgPack::SpVector::element_type& x
		, const SparseLinAlgPack::SpVector::element_type& y ) const
	{
		return ::fabs(x.value()) > ::fabs(y.value());
	}
};	// end class AbsMultVal


}	// end namespace SparseLinAlgPack 

#endif 	// SORT_BY_DESCENDING_ABS_VALUE_H
