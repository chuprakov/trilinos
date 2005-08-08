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
#include "TSFSingleScalarTypeOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "TSFRowAccessibleOp.hpp"
#include "TSFHandleable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFLinearCombination.hpp"

namespace TSFExtended
{
  using namespace TSFExtendedOps;
  /** 
   * A DiagonalOperator is a diagonal operator.
   */
  template <class Scalar> 
  class DiagonalOperator : public SingleScalarTypeOp<Scalar>,
                           public Handleable<SingleScalarTypeOp<Scalar> >,
                           public RowAccessibleOp<Scalar>
  {
  public:
    GET_RCP(SingleScalarTypeOp<Scalar>);
    /**
     * Construct a vector containing the entries on the diagonal.
     */
    DiagonalOperator(const Vector<Scalar>& vector)
      : diagonalValues_(vector),
        domain_(vector.space()),
        range_(vector.space())
    {;}


    /** Virtual dtor */
    virtual ~DiagonalOperator(){;}



    /** 
     * Apply does an element-by-element multiply between the input 
     * vector, x, and the diagonal values.
     */
    virtual void generalApply(
                       const Thyra::ETransp            M_trans
                       ,const Thyra::VectorBase<Scalar>    &x
                       ,Thyra::VectorBase<Scalar>          *y
                       ,const Scalar            alpha = 1.0
                       ,const Scalar            beta  = 0.0
                       ) const 
    {
      Vector<Scalar>  applyRes = range_.createMember();
    
      Thyra::ele_wise_prod(alpha, *(diagonalValues_.ptr().get()), x, 
                           applyRes.ptr().get());
      Thyra::Vp_StV(applyRes.ptr().get(), beta, *y);
    }



    /** Return the domain of the operator */
    RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {return domain_.ptr();}
 



    /** Return the range of the operator */
    RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > range() const 
    {return range_.ptr();}



    /** Return the kth row  */
    void getRow(const int& k, 
                Teuchos::Array<int>& indices, 
                Teuchos::Array<Scalar>& values) const
    {
      indices.resize(1);
      indices[0] = k;
      values.resize(1);
      values[0] = diagonalValues_.getElement(k);
    }



  private:

    /**
     * The vector of diagonal entries in the operator.
     */
    Vector<Scalar>  diagonalValues_;
    VectorSpace<Scalar> domain_;
    VectorSpace<Scalar> range_;
    
  };
}

#endif
