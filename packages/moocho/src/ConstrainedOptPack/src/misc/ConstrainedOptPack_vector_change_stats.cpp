// ///////////////////////////////////////////////////////
// vector_change_stats.cpp

#include <limits>

#include "../include/vector_change_stats.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

void ConstrainedOptimizationPack::vector_change_stats(
	  const VectorSlice& x, const VectorSlice& d
	, value_type* max_term, size_type* max_k
	, value_type* min_term, size_type* min_k
	, value_type* av_term )
{
	LinAlgPack::VopV_assert_sizes( x.size(), d.size() );
	const value_type
		min_num		= std::numeric_limits<value_type>::min(),
		inf			= std::numeric_limits<value_type>::infinity();
	// Initialize statistics
	*max_term	= 0.0;
	*max_k		= 0;
	*min_term	= inf;
	*min_k		= 0;
	*av_term	= 0.0;
	// Compute statistics
	VectorSlice::const_iterator
		x_itr	= x.begin(),
		d_itr	= d.begin();
	for( size_type i = 1; x_itr != x.end(); ++i, ++d_itr, ++x_itr ) {
		// Compute this ratio and make sure we are not dividing by zero.
		// We only care about ratios less than 1.0 and also deal
		// with the case that both x(i) and d(i) are zero (in which
		// case the ratio should be zero since x(i) == x(i) + d(i)).
		const value_type
			term = ::fabs(*d_itr) / ( 1 + ::fabs(*x_itr) );
		if( term > *max_term ) {
			*max_term	= term;
			*max_k		= i;
		}
		if( term < *min_term ) {
			*min_term	= term;
			*min_k		= i;
		}
		*av_term += term;
	}
	*av_term = *av_term / x.size();
}