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

#ifndef TSFSUMOPERATOR_HPP
#define TSFSUMOPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFOpDescribableByTypeID.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace TSFExtended
{
  /**
   * SumOperator is the sum (or difference) of two linear operators.
   * Whether this object represents addition or subtraction is indicated
   * by the subtraction argument to the ctor.
   */
  template <class Scalar> 
  class SumOperator : public OpDecribableByTypeID<Scalar>
  {
  public:
    /** Construct with a pair of operators and a boolean to indicate
      * if addition or substraction is to be performed. 
      */
    SumOperator(const LinearOperator& left, LinearOperator& right, bool subtraction = false)
      : left_(left.ptr()), right_(right.ptr()), subtraction_(subtraction) {;}

    /** Virtual dtor */
    virtual ~SumOperator(){;}

    /** 
     * Apply operator to a vector in the domain space and return a vector
     * in the range space.
     */
    virtual void apply(
                       const ETransp            M_trans
                       ,const TSFCore::Vector<Scalar>    &x
                       ,TSFCore::Vector<Scalar>          *y
                       ,const Scalar            alpha = 1.0
                       ,const Scalar            beta  = 0.0
                       ) const 
    {
      RefCountPtr<const TSFCore::Vector<Scalar> > applyRight;
      left_->apply(M_trans, x, y, alpha, beta);
      if (!subtraction_)
	{
	 right_->apply(M_trans, x, applyRight, alpha, 0.0); 
	}
      else 
	{
	  right_->apply(M_trans, x, applyRight, -1.0*alpha, 0.0);
	}
      linear_combination(1, 1.0, &applyRight, 1.0, y);
    }

    /** Return the domain of the operator */
    virtual RefCountPtr< const TSFCore::VectorSpace<Scalar> > domain() const {return left_->domain();}
    }

    /** Return the range of the operator */
    virtual RefCountPtr< const TSFCore::VectorSpace<Scalar> > range() const {return left_->range_();}

  private:

    RefCountPtr<const LinearOperator > left_;  
    RefCountPtr<const LinearOperator > right_; 
    bool subtraction_;

  };
}

#endif
