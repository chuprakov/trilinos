// /////////////////////////////////////////////////////////////////////
// rsqp_algo_conversion.h

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