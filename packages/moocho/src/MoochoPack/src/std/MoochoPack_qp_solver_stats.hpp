// ////////////////////////////////////////////////////////////////////
// MoochoPack_qp_solver_stats.hpp
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

#ifndef QP_SOLVER_STATS_HH
#define QP_SOLVER_STATS_HH

#include "MoochoPack_Types.hpp"
#include "ConstrainedOptPack_QPSolverStats.hpp"
#include "IterationPack_CastIQMember.hpp"

namespace MoochoPack {

/// Name given to the active set statistics iteration quantity
extern const std::string qp_solver_stats_name;

using ConstrainedOptPack::QPSolverStats;

///
/** Class for object that attempts to return an IterQuantityAccess<QPSolverStats>
  * from an AlgorithmState object with the name qp_solver_stats_name.
  */
class qp_solver_stats_iq_member
	: public CastIQMember<QPSolverStats>
{
public:
    qp_solver_stats_iq_member()
    	: CastIQMember<QPSolverStats>(qp_solver_stats_name)
    {}
};

}	// end namespace MoochoPack

#endif	// QP_SOLVER_STATS_HH
