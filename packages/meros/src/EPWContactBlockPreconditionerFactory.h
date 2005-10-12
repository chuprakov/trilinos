/*
// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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
*/

#ifndef EPWCONTACTBLOCKPRECONDITIONERFACTORY_H
#define EPWCONTACTBLOCKPRECONDITIONERFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFLinearSolver.h"
#include "TSFArray.h"
// #include "SchurFactoryBase.h"
// #include "SchurFactory.h"
#include "TSFIdentityOperator.h"
// #include "SimpleOperatorSource.h"
#include <iostream>
#include <string>

namespace TSF
{
  // using namespace TSF;
  /** \ingroup Preconditioner
   *  factory for block preconditioners for Simple
   */
  class EPWContactBlockPreconditionerFactory : public TSFPreconditionerFactoryBase
    {
    public:
      /** constructor for preconditioner using right LDU factors only */
      // Testing simpler method with two solvers on the block diagonals
      EPWContactBlockPreconditionerFactory(const TSFLinearSolver& blockSolverA,
					   const TSFLinearSolver& blockSolverC,
					   const TSFLinearSolver& blockSolverE,
					   const TSFLinearSolver& schurSolver);

      /** constructor for preconditioner using right LDU factors only */
      //EPWContactBlockPreconditionerFactory(const TSFLinearSolver& blockSolver,
      //			    const TSFLinearSolver& SchurSolver);

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& op) const;

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFOperatorSource& saddleOpSrc) const;

      /** write to a string */
      virtual string toString() const;

    private:
      TSFLinearSolver blockSolverA_;
      TSFLinearSolver blockSolverC_;
      TSFLinearSolver blockSolverE_;
      TSFLinearSolver schurSolver_;
      TSFLinearOperator op_;
    };

}
#endif
