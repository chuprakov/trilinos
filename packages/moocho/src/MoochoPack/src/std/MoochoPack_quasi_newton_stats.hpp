// ////////////////////////////////////////////////////////////////////
// quasi_newton_stats.h
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

#ifndef QUASI_NEWTON_STATS_HH
#define QUASI_NEWTON_STATS_HH

#include "QuasiNewtonStats.h"
#include "GeneralIterationPack/src/CastIQMember.h"

namespace ReducedSpaceSQPPack {

/// Name given to the quasi-Newton updating staistics iteration quantity
extern const std::string quasi_newton_stats_name;

///
/** Class for object that attempts to return an IterQuantityAccess<QuasiNewtonStats>
  * from an AlgorithmState object with the name quasi_newton_stats_name.
  */
class quasi_newton_stats_iq_member
	: public CastIQMember<QuasiNewtonStats>
{
public:
    quasi_newton_stats_iq_member()
    	: CastIQMember<QuasiNewtonStats>(quasi_newton_stats_name)
    {}
};

}	// end namespace ReducedSpaceSQPPack

#endif	// QUASI_NEWTON_STATS_HH
