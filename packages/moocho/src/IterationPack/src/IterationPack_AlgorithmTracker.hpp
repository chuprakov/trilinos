// ///////////////////////////////////////////////////////////////
// IterationPack_AlgorithmTracker.hpp
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

#include "IterationPack_Types.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace IterationPack {

///
/** Used to ouput iteration results and other information.
  *
  * This interface can be implemented by outside clients of an iterative
  * algorithm to monitor or "track" the progress of the algorithm.
  *
  * ToDo: Write more documentation!
  */
class AlgorithmTracker {
public:

  ///
  virtual ~AlgorithmTracker() {}

	/** @name Public types */
	//@{

	///
	typedef Teuchos::RefCountPtr<std::ostream>    ostream_ptr_t;

	//@}

	/** @name Constructors */
	//@{

	///
	/** Construct with an output stream for journal_out.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>journal_out.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>).
	 * </ul>
	 */
	AlgorithmTracker(const ostream_ptr_t& journal_out);
	
	//@}
	
	/** @name Algorithm iteration state notification */
	//@{
	
	///
	/** Reinitialize the track object right before it is used.
	 *
	 * The default implementation does nothing.
	 */
	virtual void initialize();

	///
	/** Output information about an iteration just completed.
	  *
	  * The default just does nothing.
	  */
	virtual void output_iteration(const Algorithm& algo) const;

	///
	/** Output information about a just completed algorithm.
	  *
	  * The default just does nothing.
	  */
	virtual void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

	//@}

	/** @name Journal file access */
	//@{

	///
	/** Set a smart pointer to the journal file.
	  */
	virtual void set_journal_out(const ostream_ptr_t& journal_out);

	///
	/** Get the smart pointer to the journal file.
	 */
	const ostream_ptr_t& get_journal_out() const;

	///
	/** Return a reference to a <tt>std::ostream</tt> to be used to output debug information 
	  * and the like.
	  */
	virtual std::ostream& journal_out() const;

	//@}

private:

#ifndef DOXYGEN_COMPILE
	ostream_ptr_t   journal_out_;
#endif

	// not defined and not to be called
	AlgorithmTracker();
	
};	// end class AlgorithmTracker

}	// end namespace IterationPack 

#endif // ALGORITHM_TRACK_H
