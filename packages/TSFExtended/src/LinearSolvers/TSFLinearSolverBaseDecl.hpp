/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFLINEARSOLVERBASEDECL_HPP
#define TSFLINEARSOLVERBASEDECL_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVector.hpp"
#include "TSFSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFLinearOperatorDecl.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /** */
  template <class Scalar>
  class LinearSolverBase 
  {
  public:
    /** */
    LinearSolverBase(const ParameterList& params);

    /** */
    virtual ~LinearSolverBase(){;}

    /** */
    virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                                      const Vector<Scalar>& rhs,
                                      Vector<Scalar>& soln) const = 0;

    /** */
    const ParameterList& parameters() const ;

    /** */
    ParameterList& parameters();

    /** */
    int getVerbosity() const ;

    /** */
    static string verbosityParam();

    /** */
    template <typename T>
    static void setParameter(const ParameterList& params,
                             T* valuePtr, 
                             const string& paramName);
  private:
    ParameterList params_;
  };
}

#endif
