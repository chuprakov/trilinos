// ////////////////////////////////////////////////////////////////////
// d_bounds_iter_quant.h
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
