// ///////////////////////////////////////////////////////////////
// AlgorithmTrack.h
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

#ifndef ALGORITHM_TRACK_H
#define ALGORITHM_TRACK_H

#include <iosfwd>

#include "GeneralIterationPackTypes.h"

namespace GeneralIterationPack {

///
/** Used to ouput iteration results and other information.
  *
  * This interface can be implemented by outside clients of an iterative
  * algorithm to monitor or "track" the progress of the algorithm.
  *
  * ToDo: Write more documentation!
  */
class AlgorithmTrack {
public:

	/** @name Constructors */
	//@{

	///
	AlgorithmTrack(std::ostream& journal_out)
		: journal_out_(&journal_out)
			{}

	//@}

	/** @name Algorithm iteration state notification */
	//@{

	///
	/** Output information about an iteration just completed.
	  *
	  * The default just does nothing.
	  */
	virtual void output_iteration(const Algorithm& algo) const
			{}

	///
	/** Output information about a just completed algorithm.
	  *
	  * The default just does nothing.
	  */
	virtual void output_final(const Algorithm& algo, EAlgoReturn algo_return) const
			{}

	//@}

	/** @name Journal file access */
	//@{

	///
	/** Set the journal file.
	  */
	virtual void set_journal_out(std::ostream& journal_out)
			{	journal_out_ = &journal_out; }

	///
	/** Return a reference to a std::ostream to be used to output debug information 
	  * and the like.
	  */
	virtual std::ostream& journal_out() const
			{	return *const_cast<AlgorithmTrack*>(this)->journal_out_; }

	//@}

protected:
	std::ostream *journal_out_;

	// not defined and not to be called
	AlgorithmTrack();
	
};	// end class AlgorithmTrack

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_TRACK_H
