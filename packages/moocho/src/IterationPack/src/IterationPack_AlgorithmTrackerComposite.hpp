// ////////////////////////////////////////////////////////////////////////////
// AlgorithmTrackComposite.h
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

#ifndef ALGORITHM_TRACK_COMPOSITE_H
#define ALGORITHM_TRACK_COMPOSITE_H

#include <list>

#include "AlgorithmTrack.h"
#include "ref_count_ptr.h"

namespace GeneralIterationPack {

///
/** This class acts a composite container for other \c AlgorithmTrack objects.
 *
 * This class exposes a <tt>std::list<AlgorithmTrack*></tt> object and lets the client
 * manipulate the list.  It is up to the client to maintain this list.
 *
 * See the "Composite" pattern in Gama et al, 1995.
 */
class AlgorithmTrackComposite : public AlgorithmTrack {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<AlgorithmTrack>      track_ptr_t;
	///
	typedef std::list<track_ptr_t>                                    track_list_t;
	///
	AlgorithmTrackComposite(std::ostream& journal_out) : AlgorithmTrack(journal_out)
	{}
	/// Give access to the list of \c AlgorithmTrack object pointers.
	track_list_t& tracks();
	///
	const track_list_t& tracks() const;
	///
	void output_iteration(const Algorithm& algo) const;
	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

private:

	track_list_t    tracks_;

};	// end class AlgorithmTrackComposite

// ///////////////////////////////////
// Inline members

inline
AlgorithmTrackComposite::track_list_t&
AlgorithmTrackComposite::tracks()
{ 
	return tracks_;
}

inline
const AlgorithmTrackComposite::track_list_t&
AlgorithmTrackComposite::tracks() const
{ 
	return tracks_;
}

}	// end namespace GeneralIterationPack 

#endif	// ALGORITHM_TRACK_COMPOSITE_H
