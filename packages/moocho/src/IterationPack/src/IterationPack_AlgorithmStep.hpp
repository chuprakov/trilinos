// ///////////////////////////////////////////////////////////////
// AlgorithmStep.h
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

#ifndef ALGORITHM_STEP_H
#define ALGORITHM_STEP_H

#include <iosfwd>

#include "GeneralIterationPackTypes.h"

namespace GeneralIterationPack {

///
/** Base type for all objects that perform steps in an <tt>Algorithm</tt>.
  */
class AlgorithmStep {
public:

	///
	typedef size_t poss_type;

	///
	/** Called by <tt>Algorithm</tt> to perform a main, pre or post step at step_poss and assoc_step_poss.
	  *
	  * @return Should return false if this step object has terminated the algorithm or
	  * redirected control to another step.  In this case it is assumed that <tt>this</tt> called
	  * <tt>algo\ref Algorithm::terminate ".terminate(...)"</tt>
	  * or <tt>algo\ref Algorithm::do_step_next ".do_step_next(...)"</tt>.
	  */
	virtual bool do_step(
		Algorithm&             algo
		,poss_type             step_poss
		,EDoStepType           type
		,poss_type             assoc_step_poss
		) = 0;

	///
	/** Called by <tt>Algorithm</tt> to inform when a runtime configuration change
	  * is finihed.
	  *
	  * The default does nothing.
	  */
	virtual void inform_updated(Algorithm& algo)
	{}

	///
	/** Called by <tt>Algorithm::print_algorithm()</tt> to print out what this step does in Matlab like format.
	  *
	  * The default does nothing.
	  */
	virtual void print_step(
		const Algorithm&         algo
		,poss_type               step_poss
		,EDoStepType             type
		,poss_type               assoc_step_poss
		, std::ostream&          out
		,const std::string&      leading_str
		) const
	{}

};	// end class AlgorithmStep

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_STEP_H
