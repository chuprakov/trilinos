// ////////////////////////////////////////////////////////////////////
// ActSetStats.h

#ifndef ACT_SET_STATS_H
#define ACT_SET_STATS_H

#include "../ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Class for storing statistics about the changes in the active set
  * of an SQP algorithm
  */
class ActSetStats {
public:

	// Public types

	/// Set to this value if a statistic is not known.
	enum { NOT_KNOWN = -1 };

	// Public interface

	/// Construct all unknowns
	ActSetStats()
		: num_active_(NOT_KNOWN), num_adds_(NOT_KNOWN), num_drops_(NOT_KNOWN)
	{}

	/// Initialize the statistics
	void set_stats( int num_active, int num_adds, int num_drops )
	{
		num_active_	= num_active;
		num_adds_	= num_adds;
		num_drops_	= num_drops;
	}

	///
	int num_active() const
	{
		return num_active_;
	}
	///
	int	num_adds() const
	{
		return num_adds_;
	}
	///
	int	num_drops() const
	{
		return num_drops_;
	}

private:
	int num_active_;
	int	num_adds_;
	int	num_drops_;

};	// end class ActSetStats

}	// end namespace ReducedSpaceSQPPack

#endif	// ACT_SET_STATS_H