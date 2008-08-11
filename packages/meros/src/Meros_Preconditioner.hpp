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

#ifndef MEROS_PRECONDITIONER_H
#define MEROS_PRECONDITIONER_H

#include "Thyra_VectorImpl.hpp" // need for LinOpDecl
#include "Thyra_VectorSpaceImpl.hpp" // need for LinOpDecl
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_RCP.hpp"


namespace Meros 
{

  using namespace Teuchos;
  using namespace Thyra;

  /** \brief 
   * Clean interface for preconditioners
   */
  template <class Scalar>
  class Preconditioner
  {
  public:
    /** */
    Preconditioner() : p_() {;}
    
    /** */
    Preconditioner(const RCP<PreconditionerBase<Scalar> >& p)
      : p_(p) {;}
    
    /** */
    Preconditioner(Handleable<PreconditionerBase<Scalar> >* p)
      : p_(p->getRcp()) {;}

    /** */
    LinearOperator<Scalar> getRight() const 
    {
      return p_->getNonconstRightPrecOp();
    }

    /** */
    bool hasRight() const 
    {
      return p_->getNonconstRightPrecOp().get() != 0;
    }

    /** */
    LinearOperator<Scalar> getLeft() const {return p_->getNonconstLeftPrecOp();}

    /** */
    bool hasLeft() const 
    {
      return p_->getNonconstLeftPrecOp().get() != 0;
    }

    /** */
    LinearOperator<Scalar> getGenericPrec() const {return p_->getNonconstGenericPrecOp();}


    /** */
    bool hasGeneric() const 
    {
      return p_->getNonconstUnspecifiedPrecOp().get() != 0;
    }

    /** */
    const RCP<PreconditionerBase<Scalar> >& ptr() const {return p_;}
    /** */
    RCP<PreconditionerBase<Scalar> > ptr() {return p_;}

  private:
    RCP<PreconditionerBase<Scalar> > p_;
  };


} // namespace Meros

#endif // MEROS_PCD_OPERATOR_SOURCE_H
