// ///////////////////////////////////////////////////////////////
// sparse_bounds.cpp

#include "../include/sparse_bounds.h"

namespace SparseLinAlgPack {

///
/** Count the number of sparse bounds where at least one bound is
  * finite.
  */
SparseLinAlgPack::size_type
SparseLinAlgPack::num_bounds( const SpVectorSlice& bl, const SpVectorSlice& bu )
{
	SpVectorSlice::const_iterator
		bl_itr			= bl.begin(),
		bl_itr_end		= bl.end(),
		bu_itr			= bu.begin(),
		bu_itr_end		= bu.end();
	size_type num_bounds = 0;
	while( bl_itr != bl_itr_end || bu_itr != bu_itr_end ) {
		if( ( bl_itr != bl_itr_end )
			&& ( bu_itr == bu_itr_end || bl_itr->indice() < bu_itr->indice() ) )
		{
			// Only the lower bound is finite
			++bl_itr;
		}
		else if( ( bu_itr != bu_itr_end )
			&& ( bl_itr == bl_itr_end || bu_itr->indice() < bl_itr->indice() ) )
		{
			// Only the upper bound is finite
			++bu_itr;
		}
		else if(bl_itr->indice() == bu_itr->indice()) {
			// Both bounds exist.
			++bl_itr;
			++bu_itr; 
		}
		++num_bounds;
	}
	return num_bounds;
}

}	// end namespace SparseLinAlgPack