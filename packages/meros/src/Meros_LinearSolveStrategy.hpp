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

#ifndef MEROS_LINEAR_SOLVE_STRATEGY_H
#define MEROS_LINEAR_SOLVE_STRATEGY_H

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorImpl.hpp" 
#include "Thyra_VectorSpaceImpl.hpp" 
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Meros_Preconditioner.hpp"

namespace Meros 
{
  using namespace Teuchos;
  using namespace Thyra;

  /** \brief 
   * Abstract specification of solve strategy
   */
  template <class Scalar>
  class LinearSolveStrategy
  {
  public:
    /** */
    LinearSolveStrategy() : lowsf_(), pf_() {;}
    
    /** */
    LinearSolveStrategy(const RCP<LinearOpWithSolveFactoryBase<Scalar> >& lowsf)
      : lowsf_(lowsf), pf_(lowsf->getPreconditionerFactory()) {;}
    
    /** */
    LinearSolveStrategy(Handleable<LinearOpWithSolveFactoryBase<Scalar> >* lowsf)
      : lowsf_(lowsf->getRcp()), pf_(lowsf_->getPreconditionerFactory()) {;}

    /** \brief . */
    RCP<LinearOpWithSolveFactoryBase<Scalar> > getLOWSF() const
      {
        return lowsf_;
      }

    /** */
    RCP<LinearOpWithSolveBase<Scalar> > getLOWS(const LinearOperator<Scalar>& op) const 
    {
      RCP<LinearOpWithSolveBase<Scalar> > rtn = lowsf_->createOp();
      RCP<const LinearOpSourceBase<Scalar> > src 
        = rcp(new DefaultLinearOpSource<Scalar>(op.ptr()));

      if (pf_.get() == 0)
        {
          LinearOpWithSolveBase<Scalar>* rtnPtr = rtn.get();
          lowsf_->initializeOp(src, rtnPtr);
        }
      else
        {
          Preconditioner<Scalar> p = pf_->createPrec();
          pf_->initializePrec(src, p.ptr().get());
          lowsf_->initializePreconditionedOp(src, p.ptr(), &*rtn);
        }

      return rtn;
    }

    /** */
    RCP<LinearOpWithSolveBase<Scalar> > getLOWS(const ConstLinearOperator<Scalar>& op) const 
    {
      RCP<LinearOpWithSolveBase<Scalar> > rtn = lowsf_->createOp();
      RCP<const LinearOpSourceBase<Scalar> > src 
        = rcp(new DefaultLinearOpSource<Scalar>(op.constPtr()));

      if (pf_.get() == 0)
        {
          LinearOpWithSolveBase<Scalar>* rtnPtr = rtn.get();
          lowsf_->initializeOp(src, rtnPtr);
        }
      else
        {
          Preconditioner<Scalar> p = pf_->createPrec();
          pf_->initializePrec(src, p.ptr().get());
          lowsf_->initializePreconditionedOp(src, p.ptr(), &*rtn);
        }

      return rtn;
    }

  protected:
    RCP<LinearOpWithSolveFactoryBase<Scalar> > lowsf() {return lowsf_;}
    RCP<PreconditionerFactoryBase<Scalar> > pf() {return pf_;}

  private:

    RCP<LinearOpWithSolveFactoryBase<Scalar> > lowsf_;
    RCP<PreconditionerFactoryBase<Scalar> > pf_;
  };


} // namespace Meros

#endif // MEROS_PCD_OPERATOR_SOURCE_H
