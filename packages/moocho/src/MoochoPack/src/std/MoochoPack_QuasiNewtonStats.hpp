// ////////////////////////////////////////////////////////////////////
// QuasiNewtonStats.h

#ifndef QUASI_NEWTON_STATS_H
#define QUASI_NEWTON_STATS_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

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

}	// end namespace ReducedSpaceSQPPack

#endif	// QUASI_NEWTON_STATS_H
