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

#ifndef MEROS_PRECONDITIONING_STRATEGY_H
#define MEROS_PRECONDITIONING_STRATEGY_H

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorImpl.hpp" // need for LinOpDecl
#include "Thyra_VectorSpaceImpl.hpp" // need for LinOpDecl
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_PreconditionerBase.hpp"

namespace Meros 
{
  using namespace Teuchos;
  using namespace Thyra;

  /** \brief 
   * Abstract specification of preconditioning strategy
   */
  template <class Scalar>
  class PreconditioningStrategy
  {
  public:
    /** */
    PreconditioningStrategy() : pf_() {;}
    
    /** */
    PreconditioningStrategy(const RCP<PreconditionerFactoryBase<Scalar> >& pf)
      : pf_(pf) {;}
    
    /** */
    PreconditioningStrategy(Handleable<PreconditionerFactoryBase<Scalar> >* pf)
      : pf_(pf->getRcp()) {;}

    /** */
    RCP<PreconditionerBase<Scalar> > getPrec(const LinearOperator<Scalar>& op) const 
    {
      RCP<PreconditionerBase<Scalar> > rtn = pf_->createPrec();
      pf_->initializeOp<Scalar>(op.ptr(), &*rtn);
      return rtn;
    }

    /** */
    const RCP<PreconditionerFactoryBase<Scalar> >& ptr() const {return pf_;}
    /** */
    RCP<PreconditionerFactoryBase<Scalar> > ptr() {return pf_;}

  protected:
    RCP<PreconditionerFactoryBase<Scalar> > pf() {return pf_;}

  private:

    RCP<PreconditionerFactoryBase<Scalar> > pf_;
  };


} // namespace Meros

#endif // MEROS_PCD_OPERATOR_SOURCE_H
