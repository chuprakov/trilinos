// ///////////////////////////////////////////////////////////////////////
// ReducedSpaceSQPPackExceptions.h

#ifndef REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
#define REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H

#include "ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

// Thrown if the constraints are infeasible
class InfeasibleConstraints : public std::logic_error
{public: InfeasibleConstraints(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a line search failure occurs.
class LineSearchFailure : public std::logic_error
{public: LineSearchFailure(const std::string& what_arg) : std::logic_error(what_arg){}};

}	// end namespace ReducedSpaceSQPPack 

#endif // REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H