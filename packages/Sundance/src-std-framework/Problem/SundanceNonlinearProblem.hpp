/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_NONLINEARPROBLEM_H
#define SUNDANCE_NONLINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFNonlinearOperatorBase.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

    /** 
     * NonlinearProblem encapsulates a discrete nonlinear problem, and can
     * be passed to a nonlinear solver such as NOX.
     */
  class NonlinearProblem 
    : public TSFExtended::ObjectWithVerbosity<NonlinearProblem>,
      public TSFExtended::NonlinearOperatorBase<double>
    {
    public:
      /** Empty ctor */
      NonlinearProblem();

      /** Construct with a mesh, equation set, bcs, test and unknown funcs,
       * and a vector type */
      NonlinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
                       const Expr& test, const Expr& unk, const Expr& u0, 
        const TSFExtended::VectorType<double>& vecType,
        bool partitionBCs = false);

      /** Construct with a mesh, equation set, bcs, test and unknown funcs,
       * parameters, and a vector type */
      NonlinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
                       const Expr& test, const Expr& unk, const Expr& u0, 
                       const Expr& params, const Expr& paramVals,  
        const TSFExtended::VectorType<double>& vecType,
        bool partitionBCs = false);


#ifndef DOXYGEN_DEVELOPER_ONLY
      /** */
      NonlinearProblem(const RefCountPtr<Assembler>& assembler, 
                       const Expr& u0);
#endif

      /** Compute the residual and Jacobian at the current evaluation point */
      LinearOperator<double> computeJacobianAndFunction(Vector<double>& functionValue) const ;

      /** Write the Jacobian and residual into the objects provided */
      void computeJacobianAndFunction(LinearOperator<double>& J,
                                      Vector<double>& resid) const ;

      /** Return the current evaluation point as a Sundance expression */
      Expr getU0() const {return u0_;}
      
      /** Compute the residual at the current eval point */
      TSFExtended::Vector<double> computeFunctionValue() const ;
      
      /** Write the residual into the object provided */
      void computeFunctionValue(Vector<double>& resid) const ;

      /** Get an initial guess */
      TSFExtended::Vector<double> getInitialGuess() const ;

      /** Create the Jacobian object, but don't fill it in. */
      LinearOperator<double> allocateJacobian() const ;
      



#ifndef DOXYGEN_DEVELOPER_ONLY

      /* Handle boilerplate */
      GET_RCP(TSFExtended::NonlinearOperatorBase<double>);

    private:
      
      /** */
      RefCountPtr<Assembler> assembler_;

      /** */
      mutable TSFExtended::LinearOperator<double> J_;

      /** */
      Expr u0_;

      /** */
      mutable DiscreteFunction* discreteU0_;
#endif /* DOXYGEN_DEVELOPER_ONLY */
      
    };
}




#endif
