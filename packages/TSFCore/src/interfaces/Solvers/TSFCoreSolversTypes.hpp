// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////////////////
// TSFCoreSolversTypes.hpp

#ifndef TSF_CORE_SOLVERS_TYPES_HPP
#define TSF_CORE_SOLVERS_TYPES_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {
namespace Solvers {

///
const int DEFAULT_MAX_ITER = -1;

///
enum ESolveReturn {
	SOLVED_TO_TOL          ///< The linear system(s) solved to tolerance
	,MAX_ITER_EXCEEDED     ///< <tt>max_iter</tt> exceeded
};

///
struct SolveReturn {
	SolveReturn(ESolveReturn solve_status_in, int num_iter_in) :solve_status(solve_status_in),num_iter(num_iter_in) {} 
	ESolveReturn    solve_status;
	int             num_iter;
};

namespace Exceptions {

/** \defgroup TSFCoreSolversExceptions_grp Basic TSFCore::Solvers exception types.
 */
//@{

/// Thrown if a linear solve failed
class FailureToConverge : public std::logic_error
{public: FailureToConverge(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a the operator is defective in some way (related to linear solver method)
class SolverBreakdown : public FailureToConverge
{public: SolverBreakdown(const std::string& what_arg) : FailureToConverge(what_arg) {}};

/// Thrown if operator turns out to be indefinite
class Indefinite : public SolverBreakdown
{public: Indefinite(const std::string& what_arg) : SolverBreakdown(what_arg) {}};

//@}

} // namespace Exceptions

template<class Scalar> class SolverState;
template<class Scalar> class ConvergenceTester;
template<class Scalar> class IterativeLinearSolver;

} // namespace Solvers
} // namespace TSFCore

#endif // TSF_CORE_SOLVERS_TYPES_HPP
