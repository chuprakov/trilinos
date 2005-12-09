// ///////////////////////////////////////////////////////
// vector_change_stats.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <limits>

#include "ConstrainedOptPack_vector_change_stats.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

void ConstrainedOptPack::vector_change_stats(
	  const DVectorSlice& x, const DVectorSlice& d
	, value_type* max_term, size_type* max_k
	, value_type* min_term, size_type* min_k
	, value_type* av_term )
{
	DenseLinAlgPack::VopV_assert_sizes( x.size(), d.size() );
	const value_type
		min_num		= std::numeric_limits<value_type>::min(),
		inf			= std::numeric_limits<value_type>::max();
	// Initialize statistics
	*max_term	= 0.0;
	*max_k		= 0;
	*min_term	= inf;
	*min_k		= 0;
	*av_term	= 0.0;
	// Compute statistics
	DVectorSlice::const_iterator
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
