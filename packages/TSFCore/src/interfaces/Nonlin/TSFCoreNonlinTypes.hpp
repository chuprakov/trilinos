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

// ///////////////////////////////////////////////////////////////
// TSFCoreNonlinTypes.hpp

#ifndef TSF_CORE_NONLIN_TYPES_HPP
#define TSF_CORE_NONLIN_TYPES_HPP

#include "TSFCoreSolversTypes.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Enumeration for determining if an operator is nonsingular or not.
 */
enum ENonsingStatus {
	OP_NONSINGULAR            ///< The operator is known to be (or should be) nonsingular
	,OP_SINGULAR              ///< The oeprator is known to be (or should be) singular
	,OP_SINGULARITY_UNKNOWN   ///< Unknown if the operator is singular or nonsingular
};

template<class Scalar> class LinearSolveOp;
template<class Scalar> class LinearOpWithSolve;
template<class Scalar> class NonlinearProblem;
template<class Scalar> class NonlinearProblemFirstOrder;

} // namespace Nonlin
} // namespace TSFCore

#endif // TSF_CORE_NONLIN_TYPES_HPP
