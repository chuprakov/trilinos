// ////////////////////////////////////////////////////////////////////////
// AlgorithmTrackerComposite.cpp
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

#include <algorithm>

#include "IterationPack_AlgorithmTrackerComposite.hpp"

namespace IterationPack {

AlgorithmTrackerComposite::AlgorithmTrackerComposite(const ostream_ptr_t& journal_out)
	: AlgorithmTracker(journal_out)
{}

void AlgorithmTrackerComposite::initialize()
{
	track_list_t::const_iterator
		itr = tracks_.begin(), itr_end = tracks_.end();
	for(; itr != itr_end; ++itr)
		(*itr)->initialize();
}

void AlgorithmTrackerComposite::output_iteration(
	const Algorithm& algo
	) const
{
	track_list_t::const_iterator
		itr = tracks_.begin(), itr_end = tracks_.end();
	for(; itr != itr_end; ++itr)
		(*itr)->output_iteration(algo);
}

void AlgorithmTrackerComposite::output_final(
	const Algorithm& algo, EAlgoReturn algo_return
	) const
{
	track_list_t::const_iterator
		itr = tracks_.begin(), itr_end = tracks_.end();
	for(; itr != itr_end; ++itr)
		(*itr)->output_final(algo,algo_return);
}

} // end namespace IterationPack
