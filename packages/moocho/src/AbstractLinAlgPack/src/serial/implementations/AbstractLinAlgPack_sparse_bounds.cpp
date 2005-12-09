// ///////////////////////////////////////////////////////////////
// sparse_bounds.cpp
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

#include "AbstractLinAlgPack_sparse_bounds.hpp"

///
/** Count the number of sparse bounds where at least one bound is
  * finite.
  */
AbstractLinAlgPack::size_type
AbstractLinAlgPack::num_bounds( const SpVectorSlice& bl, const SpVectorSlice& bu )
{
	SpVectorSlice::const_iterator
		bl_itr			= bl.begin(),
		bl_itr_end		= bl.end(),
		bu_itr			= bu.begin(),
		bu_itr_end		= bu.end();
	size_type num_bounds = 0;
	while( bl_itr != bl_itr_end || bu_itr != bu_itr_end ) {
		if( ( bl_itr != bl_itr_end )
			&& ( bu_itr == bu_itr_end || bl_itr->indice() + bl.offset() < bu_itr->indice() + bu.offset() ) )
		{
			// Only the lower bound is finite
			++bl_itr;
		}
		else if( ( bu_itr != bu_itr_end )
			&& ( bl_itr == bl_itr_end || bu_itr->indice() + bu.offset() < bl_itr->indice() + bl.offset()) )
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
