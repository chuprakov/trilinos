// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_MoochoTrackerStatsStd.hpp
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

#ifndef RSQP_TRACK_STATS_STD_H
#define RSQP_TRACK_STATS_STD_H

#include "MoochoPack_quasi_newton_stats.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "StopWatchPack_stopwatch.hpp"

namespace MoochoPack {

///
/** This is a simple track class for getting statistics about a solved (or not
 * solved) NLP.
 *
 * When output_final(...) is called the output stream will have the following
 * written to it:
 *
 \begin{verbatim}
 status         =      solved; # solved, except, max_iter, max_run_time
 niter          =          10; # Number of rSQP iterations
 nfunc          =          15; # max( number f(x) evals, number c(x) evals ) 
 ngrad          =          13; # max( number Gf(x) evals, number Gc(x) evals ) 
 CPU            =        0.50; # Number of CPU seconds total
 obj_func       =    1.046e-2; # Objective function value f(x) at final point
 feas_kkt_err   =   2.457e-10; # Feasibility error at final point (scaled ||c(x)||inf)
 opt_kkt_err    =    4.568e-7; # Optimality error at final point (scaled ||rGL||inf)
 nact           =          40; # Number of total active constraints at the final point
 nbasis_change  =           1; # Number of basis changes
 nquasi_newton  =           6; # Number of quasi-newton updates
 \end{verbatim}
 *
 * Any statistic that is not known will be given the value '-'.  If the returned status
 * is 'execpt' then some exception was thrown or some other error occured so current
 * information may not be available.  In this case every effort is made to fill the rest
 * of the information from prior iterations.
 * The names of these fields will not change and the 'stat = value; # comment' format
 * can be counted.  However, the spacing and the precision of the numbers may be different
 * from what is shown above.
 */
class MoochoTrackerStatsStd
	: public IterationPack::AlgorithmTracker
 {
public:

	/// Construct with an output stream object.
	MoochoTrackerStatsStd( const ostream_ptr_t& o, const ostream_ptr_t& journal_out );

	///
	/* Set the output stream for statistics outputting.
	 */
	void set_output_stream(const ostream_ptr_t& o);

	/// Get the output stream for statistics outputting.
	const ostream_ptr_t& get_output_stream() const;

	/** @name Overridden from AlgorithmTracker */
	//@{

	/// Restarts the timer
	void initialize();
	///
	void output_iteration(const Algorithm& algo) const;
	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

	//@}

protected:

	std::ostream& o() const
	{	return *o_; }

private:
	ostream_ptr_t                       o_;
	mutable int		                    num_QN_updates_;
	quasi_newton_stats_iq_member	    quasi_newton_stats_;
	mutable StopWatchPack::stopwatch    timer_;

	// Not defined and not to be called
	MoochoTrackerStatsStd();

};	// end class MoochoTrackerStatsStd

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_STATS_STD_H
