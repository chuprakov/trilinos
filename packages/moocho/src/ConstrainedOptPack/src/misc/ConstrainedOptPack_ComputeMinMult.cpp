// ////////////////////////////////////////////////////////
// ComputeMinMult.cpp
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

#include "ConstrainedOptPack_ComputeMinMult.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

namespace {
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace

ConstrainedOptPack::value_type
ConstrainedOptPack ::min_abs( const DVectorSlice& mu )
{
	if( !mu.dim() )
		return 0.0;
	value_type min = ::fabs(mu(1));
	for( DVectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); )
		min = my_min( min, ::fabs(*itr++) );
	return min;
}

ConstrainedOptPack::value_type
ConstrainedOptPack ::min_abs( const SpVectorSlice& mu )
{
	if( !mu.dim() )
		return 0.0;
	if( !mu.nz() )
		return 0.0;
	value_type min = ::fabs(mu.begin()->value());
	for( SpVectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); ++itr )
		min = my_min( min, ::fabs(itr->value()) );
	return min;
}
