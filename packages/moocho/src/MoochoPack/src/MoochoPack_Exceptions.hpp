// ///////////////////////////////////////////////////////////////////////
// MoochoPack_Exceptions.hpp
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

#ifndef REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
#define REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H

#include "MoochoPack_Types.hpp"
#include "ConstrainedOptPack_QPSolverStats.hpp"

namespace MoochoPack {

/** \defgroup MoochoPack_grp Standard exceptions for MoochoPack */
//@{

// Thrown if the constraints are infeasible
class InfeasibleConstraints : public std::logic_error
{public: InfeasibleConstraints(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a line search failure occurs.
class LineSearchFailure : public std::runtime_error
{public: LineSearchFailure(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a runtime test failed.
class TestFailed : public std::runtime_error
{public: TestFailed(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a the QP failed and was not corrected
class QPFailure : public std::runtime_error
{
public:
	QPFailure(const std::string& what_arg
			  , const ConstrainedOptPack::QPSolverStats& _qp_stats)
		: std::runtime_error(what_arg)
		, qp_stats(_qp_stats)
		{}
	ConstrainedOptPack::QPSolverStats qp_stats;
};

//@}

}	// end namespace MoochoPack 

#endif // REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
