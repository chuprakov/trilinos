// ////////////////////////////////////////////////////////////////////
// MoochoPack_QuasiNewtonStats.hpp
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

#ifndef QUASI_NEWTON_STATS_H
#define QUASI_NEWTON_STATS_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

///
/** Class for storing statistics about the Quasi-Newton updating
  */
class QuasiNewtonStats {
public:

	// Public types

	/// Set to this value if a statistic is not known.
	enum EUpdate { UNKNOWN, REINITIALIZED, UPDATED, DAMPENED_UPDATED
		, SKIPED, INDEF_SKIPED };

	// Public interface

	/// Construct all unknowns
	QuasiNewtonStats()
		: update_(UNKNOWN)
	{}

	/// Initialize the statistics
	void set_updated_stats( EUpdate update )
	{
		update_ = update;
	}

	///
	EUpdate updated() const
	{
		return update_;
	}

private:
	EUpdate update_;

};	// end class QuasiNewtonStats

}	// end namespace MoochoPack

#endif	// QUASI_NEWTON_STATS_H
