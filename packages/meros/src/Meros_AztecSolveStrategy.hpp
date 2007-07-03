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

#ifndef MEROS_AZTEC_SOLVE_STRATEGY_H
#define MEROS_AZTEC_SOLVE_STRATEGY_H

#include "Meros_LinearSolveStrategy.hpp"
#include "Meros_PreconditioningStrategy.hpp"
#include "Meros_Preconditioner.hpp"
#include "AztecOO.h"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"


namespace Meros 
{
  using namespace Teuchos;
  using namespace Thyra;

  /** \brief 
   * Specification of an Aztec solve strategy
   */
  class AztecSolveStrategy
    : public AztecOOLinearOpWithSolveFactory,
      public Handleable<LinearOpWithSolveFactoryBase<double> >
  {
  public:
    /** */
    AztecSolveStrategy(const ParameterList& params) 
      : AztecOOLinearOpWithSolveFactory()
    {
      RCP<ParameterList> p = rcp(new ParameterList(params));
      this->setParameterList(p);
    }

    /** */
    AztecSolveStrategy(const ParameterList& params, 
                       const PreconditioningStrategy<double>& prec) 
      : AztecOOLinearOpWithSolveFactory()
    {
      RCP<ParameterList> p = rcp(new ParameterList(params));
      this->setParameterList(p);
      this->setPreconditionerFactory(prec.ptr(), "dummy name");
    }

    TEUCHOS_GET_RCP(LinearOpWithSolveFactoryBase<double>);

    /** */
    static ParameterList getParameters()
    {
      AztecOOLinearOpWithSolveFactory dummy;
      return *dummy.getValidParameters();
    } 

    

  };


} // namespace Meros

#endif // MEROS_PCD_OPERATOR_SOURCE_H
