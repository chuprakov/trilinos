// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_MoochoTrackerSummaryStd.hpp
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

#ifndef RSQP_TRACK_SUMMARY_STD_H
#define RSQP_TRACK_SUMMARY_STD_H

#include "MoochoPack_quasi_newton_stats.hpp"
#include "MoochoPack_qp_solver_stats.hpp"
#include "MoochoPack_act_set_stats.hpp"
#include "IterationPack_AlgorithmTracker.hpp"

namespace MoochoPack {

///
/** This class simply outputs the convergence information
  * for each iteration.
  */
class MoochoTrackerSummaryStd
	: public IterationPack::AlgorithmTracker
{
public:

	///
	enum EOptError { OPT_ERROR_REDUCED_GRADIENT_LAGR, OPT_ERROR_GRADIENT_LAGR };

	/// Construct with an output stream
	MoochoTrackerSummaryStd(
		const ostream_ptr_t      &o
		,const ostream_ptr_t     &journal_out
		,EOptError               opt_error = OPT_ERROR_REDUCED_GRADIENT_LAGR
		);

	/// Set the output stream for summary outputting
	void set_output_stream(const ostream_ptr_t& o);

	/// Get the output stream for summary outputting.
	const ostream_ptr_t& get_output_stream() const;

	///
	/** Output the total number of qp iterations back to and
	  * the k=0 iteration.
	  */
	int num_total_qp_iter() const
	{	return num_total_qp_iter_;	}

	/** @name Overridden from AlgorithmTracker */
	//@{

	///
	void output_iteration(const Algorithm& algo) const;
	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;
	
	//@}

protected:

	/// Print the header to the output
	void print_header(const NLPAlgoState &s) const;

	std::ostream& o() const
	{	return *o_; }

private:
	ostream_ptr_t                  o_;
	EOptError                      opt_error_;
	mutable int                     num_total_qp_iter_;
	quasi_newton_stats_iq_member    quasi_newton_stats_;
	qp_solver_stats_iq_member       qp_solver_stats_;
	act_set_stats_iq_member		    act_set_stats_;

	// Not defined and not to be called
	MoochoTrackerSummaryStd();

};	// end class MoochoTrackerSummaryStd

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_SUMMARY_STD_H
