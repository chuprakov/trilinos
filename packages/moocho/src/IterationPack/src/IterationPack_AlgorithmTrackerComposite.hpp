// ////////////////////////////////////////////////////////////////////////////
// IterationPack_AlgorithmTrackerComposite.hpp
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

#include "IterationPack_AlgorithmTracker.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace IterationPack {

///
/** This class acts a composite container for other \c AlgorithmTracker objects.
 *
 * This class exposes a <tt>std::list<AlgorithmTracker*></tt> object and lets the client
 * manipulate the list.  It is up to the client to maintain this list.
 *
 * See the "Composite" pattern in "Design Patterns", Gama et al, 1995.
 */
class AlgorithmTrackerComposite : public AlgorithmTracker {
public:

	///
	typedef Teuchos::RefCountPtr<AlgorithmTracker>      track_ptr_t;
	///
	typedef std::list<track_ptr_t>                                    track_list_t;
	///
	AlgorithmTrackerComposite(const ostream_ptr_t& journal_out);
	/// Give access to the list of \c AlgorithmTracker object pointers.
	track_list_t& tracks();
	///
	const track_list_t& tracks() const;

	/**  @name Overridden from AlgorithmTracker */
	//@{

	///
	void initialize();
	///
	void output_iteration(const Algorithm& algo) const;
	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	AlgorithmTracker  *tracks;
#else
	track_list_t    tracks_;
#endif

};	// end class AlgorithmTrackerComposite

// ///////////////////////////////////
// Inline members

inline
AlgorithmTrackerComposite::track_list_t&
AlgorithmTrackerComposite::tracks()
{ 
	return tracks_;
}

inline
const AlgorithmTrackerComposite::track_list_t&
AlgorithmTrackerComposite::tracks() const
{ 
	return tracks_;
}

}	// end namespace IterationPack 

#endif	// ALGORITHM_TRACK_COMPOSITE_H
