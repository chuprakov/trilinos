// ////////////////////////////////////////////////////////////////////
// MoochoPack_ActSetStats.hpp
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

#ifndef ACT_SET_STATS_H
#define ACT_SET_STATS_H

#include "../MoochoPack_Types.hpp"

namespace MoochoPack {

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
		, num_active_indep_(NOT_KNOWN), num_adds_indep_(NOT_KNOWN), num_drops_indep_(NOT_KNOWN)
	{}

	/// Initialize the statistics
	void set_stats(
		int num_active, int num_adds, int num_drops
		,int num_active_indep, int num_adds_indep, int num_drops_indep
		)
	{
		num_active_        = num_active;
		num_adds_          = num_adds;
		num_drops_         = num_drops;
		num_active_indep_  = num_active_indep;
		num_adds_indep_    = num_adds_indep;
		num_drops_indep_   = num_drops_indep;
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

	///
	int num_active_indep() const
	{
		return num_active_indep_;
	}
	///
	int	num_adds_indep() const
	{
		return num_adds_indep_;
	}
	///
	int	num_drops_indep() const
	{
		return num_drops_indep_;
	}

private:
	int num_active_;
	int	num_adds_;
	int	num_drops_;
	int num_active_indep_;
	int	num_adds_indep_;
	int	num_drops_indep_;

};	// end class ActSetStats

}	// end namespace MoochoPack

#endif	// ACT_SET_STATS_H
