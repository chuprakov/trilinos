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

#ifndef TSFDIAGONALOPERATOR_HPP
#define TSFDIAGONALOPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFOpDescribableByTypeID.hpp"
#include "TSFHandleable.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace TSFExtended
{
  /** 
   * A DiagonalOperator is a diagonal operator.
   */
  template <class Scalar> 
  class DiagonalOperator : public OpDescribableByTypeID<Scalar>,
			   public Handleable<TSFCore::LinearOp<Scalar> >
  {
  public:
    GET_RCP(TSFCore::LinearOp<Scalar>);
    /**
     * Construct a vector containing the entries on the diagonal.
     */
    DiagonalOperator(const TSFCore::Vector<Scalar>& diagonalValues)
      : diagonalValues_(diagonalValues.ptr()) {;}

    /** Virtual dtor */
    virtual ~DiagonalOperator(){;}

    /** 
     * Apply does an element-by-element multiply between the input 
     * vector, x, and the diagonal values.
     */
    virtual void apply(
                       const TSFCore::ETransp            M_trans
                       ,const TSFCore::Vector<Scalar>    &x
                       ,TSFCore::Vector<Scalar>          *y
                       ,const Scalar            alpha = 1.0
                       ,const Scalar            beta  = 0.0
                       ) const 
    {
      RefCountPtr<TSFCore::Vector<Scalar> > applyRes;
      ele_wise_product(alpha, *diagonalValues_, x, applyRes);
      Vp_StV(applyRes, beta, *y);
    }

    /** Return the domain of the operator */
    virtual RefCountPtr< const TSFCore::VectorSpace<Scalar> > domain() const {return diagonalValues_->space();}
 

    /** Return the range of the operator */
    virtual RefCountPtr< const TSFCore::VectorSpace<Scalar> > range() const {return diagonalValues_->space();}

  private:

    /**
     * The vector of diagonal entries in the operator.
     */
    RefCountPtr<const TSFCore::Vector<Scalar> > diagonalValues_;
    
  };
}

#endif
