// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_MoochoTrackerXMLSummary.hpp
//
// Copyright (C) 2001
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

#ifndef MOOCHO_TRACKER_XML_SUMMARY_HPP
#define MOOCHO_TRACKER_XML_SUMMARY_HPP

#include "IterationPack_Algorithm.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "MoochoPack_Types.hpp"

namespace MoochoPack {


using IterationPack::Algorithm;
using IterationPack::EAlgoReturn;

///
/** This class outputs an XML summary file of the algorithm 
 *   results and performance
 */
class MoochoTrackerXMLSummary
	: public IterationPack::AlgorithmTracker
{
public:

	/// Construct with an output stream
	MoochoTrackerXMLSummary(
	  const Teuchos::RefCountPtr<std::ostream> &journal_out
	  ,const std::string xml_filename
	  ,const std::string problem_name
	  ,const std::string algorithm_description
	  );

	/// Set the output stream for summary outputting
	//void set_output_stream(const ostream_ptr_t& o);

	/// Get the output stream for summary outputting.
	//const ostream_ptr_t& get_output_stream() const;

	/// Output a basic file (with failed status)
	//   that will be overwritten if there is no
	//   exception
	void output_pre_file() const;

	/** @name Overridden from AlgorithmTracker */
	//@{

	///
	void output_iteration(const Algorithm& algo) const;
	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;
	
	//@}

protected:

	/// Print the header to the output
	void open_problem_element( std::ostream& out, const Algorithm& algo) const;
	void close_problem_element( std::ostream& out) const;

private:

	mutable value_type obj_value_;
	mutable value_type c_norm_value_;

	std::string xml_filename_;
	std::string problem_name_;
	std::string algorithm_description_;

	// Not defined and not to be called
	MoochoTrackerXMLSummary();

};	// end class MoochoTrackerXMLSummary

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_SUMMARY_STD_H
