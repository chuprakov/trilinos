// ////////////////////////////////////////////////////////////////////
// MoochoPack_act_set_stats.hpp
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

#ifndef ACT_SET_STATS_HH
#define ACT_SET_STATS_HH

#include "MoochoPack_ActSetStats.hpp"
#include "IterationPack_CastIQMember.hpp"

namespace MoochoPack {

/// Name given to the active set statistics iteration quantity
extern const std::string act_set_stats_name;

///
/** Class for object that attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name act_set_stats_name.
  */
class act_set_stats_iq_member
	: public CastIQMember<ActSetStats>
{
public:
    act_set_stats_iq_member()
    	: CastIQMember<ActSetStats>(act_set_stats_name)
    {}
};

}	// end namespace MoochoPack

#endif	// ACT_SET_STATS_HH
