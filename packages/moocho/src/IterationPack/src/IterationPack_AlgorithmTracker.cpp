// ///////////////////////////////////////////////////////////////
// AlgorithmTrack.cpp
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

#include "GeneralIterationPack/src/AlgorithmTrack.h"

namespace GeneralIterationPack {

AlgorithmTrack::AlgorithmTrack(const ostream_ptr_t& journal_out)
	: journal_out_(journal_out)
{}

void AlgorithmTrack::initialize()
{}
	
void AlgorithmTrack::output_iteration(const Algorithm& algo) const
{}

void AlgorithmTrack::output_final(const Algorithm& algo, EAlgoReturn algo_return) const
{}

void AlgorithmTrack::set_journal_out(const ostream_ptr_t& journal_out)
{
	journal_out_ = journal_out;
}

const AlgorithmTrack::ostream_ptr_t&
AlgorithmTrack::get_journal_out() const
{
	return journal_out_;
}

std::ostream&
AlgorithmTrack::journal_out() const
{
	return *journal_out_;
}

}	// end namespace GeneralIterationPack 
