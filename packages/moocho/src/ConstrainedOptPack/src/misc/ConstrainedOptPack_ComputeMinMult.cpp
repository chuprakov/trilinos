// ////////////////////////////////////////////////////////
// ComputeMinMult.cpp

#include "../include/ComputeMinMult.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/VectorClass.h"

ConstrainedOptimizationPack::value_type
ConstrainedOptimizationPack ::min_abs( const VectorSlice& mu )
{
	if( !mu.size() )
		return 0.0;
	value_type min = ::fabs(mu(1));
	for( VectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); )
		min = std::_MIN( min, ::fabs(*itr++) );
	return min;
}

ConstrainedOptimizationPack::value_type
ConstrainedOptimizationPack ::min_abs( const SpVectorSlice& mu )
{
	if( !mu.size() )
		return 0.0;
	if( !mu.nz() )
		return 0.0;
	value_type min = ::fabs(mu.begin()->value());
	for( SpVectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); ++itr )
		min = std::_MIN( min, ::fabs(itr->value()) );
	return min;
}
