// /////////////////////////////////////////////////////////////////////
// rsqp_algo_conversion.h
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

#include "rSQPAlgo.h"
#include "GeneralIterationPack/include/Algorithm.h"

namespace ReducedSpaceSQPPack {

/// Convert from a Algorithm to a rSQPAlgo
inline
rSQPAlgo& rsqp_algo(Algorithm& algo)
{	return dynamic_cast<rSQPAlgo&>(algo); }

/// Convert from a Algorithm to a rSQPAlgo
inline
const rSQPAlgo& rsqp_algo(const Algorithm& algo)
{	return dynamic_cast<const rSQPAlgo&>(algo); }

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_CONVERSION_H
