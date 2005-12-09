// /////////////////////////////////////////////////////////////////////
// MoochoPack_moocho_algo_conversion.hpp
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

#ifndef RSQP_ALGO_CONVERSION_H
#define RSQP_ALGO_CONVERSION_H

#include "MoochoPack_NLPAlgo.hpp"
#include "IterationPack_Algorithm.hpp"

namespace MoochoPack {

/// Convert from a Algorithm to a NLPAlgo
inline
NLPAlgo& rsqp_algo(Algorithm& algo)
{	return dynamic_cast<NLPAlgo&>(algo); }

/// Convert from a Algorithm to a NLPAlgo
inline
const NLPAlgo& rsqp_algo(const Algorithm& algo)
{	return dynamic_cast<const NLPAlgo&>(algo); }

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_CONVERSION_H
