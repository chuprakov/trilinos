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

#include "ConstrainedOptimizationPack/include/ComputeMinMult.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/VectorClass.h"

ConstrainedOptimizationPack::value_type
ConstrainedOptimizationPack ::min_abs( const VectorSlice& mu )
{
	if( !mu.dim() )
		return 0.0;
	value_type min = ::fabs(mu(1));
	for( VectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); )
		min = std::_MIN( min, ::fabs(*itr++) );
	return min;
}

ConstrainedOptimizationPack::value_type
ConstrainedOptimizationPack ::min_abs( const SpVectorSlice& mu )
{
	if( !mu.dim() )
		return 0.0;
	if( !mu.nz() )
		return 0.0;
	value_type min = ::fabs(mu.begin()->value());
	for( SpVectorSlice::const_iterator itr = mu.begin() + 1; itr != mu.end(); ++itr )
		min = std::_MIN( min, ::fabs(itr->value()) );
	return min;
}
