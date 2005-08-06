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
#include "TSFSingleScalarTypeOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "TSFOpDescribableByTypeID.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFHandleable.hpp"
#include "TSFRowAccessibleOp.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"



namespace TSFExtended
{
  /**
   * SumOperator is the sum (or difference) of two linear operators.
   * Whether this object represents addition or subtraction is indicated
   * by the subtraction argument to the ctor.
   */
  template <class Scalar> 
  class SumOperator : public SingleScalarTypeOp<Scalar>,
                      public Handleable<SingleScalarTypeOp<Scalar> >,
                      public RowAccessibleOp<Scalar>
  {
  public:
    GET_RCP(SingleScalarTypeOp<Scalar>);

    /** Construct with a pair of operators and a boolean to indicate
     * if addition or substraction is to be performed. 
     */
    SumOperator(const LinearOperator<Scalar>& left, 
                const LinearOperator<Scalar>& right, 
                bool subtraction = false)
      : left_(left.ptr()), 
        right_(right.ptr()), 
        subtraction_(subtraction),
        domain_(left.domain()),
        range_(left.range())
    {
      TEST_FOR_EXCEPTION(domain_ != right_.domain() || range_ != right.range(),
                         runtime_error,
                         "SumOperator cannot be formed since the range "
                         << "and domain spaces of the two operators are "
                         << "not the same");
    }

    /** Virtual dtor */
    virtual ~SumOperator(){;}

    /** 
     * Apply operator to a vector in the domain space and return a vector
     * in the range space.
     */
    virtual void generalApply(
                       Thyra::ETransp            M_trans
                       ,const Thyra::VectorBase<Scalar>    &x
                       ,Thyra::VectorBase<Scalar>          *y
                       ,const Scalar            alpha = 1.0
                       ,const Scalar            beta  = 0.0
                       ) const 
    {
      Vector<Scalar>  applyRight = range_.createMember();
      left_.ptr()->generalApply(M_trans, x, y, alpha, beta);
      if (!subtraction_)
        {
          right_.ptr()->generalApply(M_trans, x, applyRight.ptr().get(), alpha, 0.0); 
        }
      else 
        {
          right_.ptr()->generalApply(M_trans, x, applyRight.ptr().get(), 
                              -1.0*alpha, 0.0);
        }
      //Thyra::linear_combination(1, 1.0, &applyRight, 1.0, y);
      Thyra::Vp_StV(y, 1.0, *(applyRight.ptr().get()));
    }

    /** Return the domain of the operator */
    virtual RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {return domain_.ptr();}
  

    /** Return the range of the operator */
    virtual RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > range() const 
    {return range_.ptr();}


    /** Return the kth row  */
    void getRow(const int& row, 
                Teuchos::Array<int>& indices, 
                Teuchos::Array<Scalar>& values) const
    {
      /* Get the row for left_ and the row for right_ and combing */
      indices.resize(0);
      values.resize(0);
      Teuchos::Array<int> indL;
      Teuchos::Array<int> indR;

      Teuchos::Array<double> valL;
      Teuchos::Array<double> valR;

      left_.getRow(row, indL, valL);
      right_.getRow(row, indR, valR);

      sort(indL, valL);
      sort(indR, valR);

      int kR = 0;
      int kL = 0;
      int L = 0;
      int R = 0;
      double vL = 0.0;
      double vR = 0.0;
      getNext(kR, R, vR, indR, valR);
      getNext(kL, L, vL, indL, valL);

      int large = indL.size() + indR.size();
      /*  Main loop to combing the arrays  */
      for (int i = 0; i < large; i++)
        {
          if (L == -1 && R == -1) break;
          else if (L == -1 || (R != -1 && L > R))
            {
              indices.append(R);
              values.append(vR);
              getNext(kR, R, vR, indR, valR);
            }
          else if (R == -1 || L < R)
            {
              indices.append(L);
              values.append(vL);
              getNext(kL, L, vL, indL, valL);
            }
          else if (L == R)
            {
              indices.append(L);
              if (subtraction_)
                {
                  values.append(vL - vR);
                }
              else
                {
                  values.append(vL + vR);
                }
              getNext(kR, R, vR, indR, valR);
              getNext(kL, L, vL, indL, valL);
            }
        }
    }



  private:
    LinearOperator<Scalar>  left_;  
    LinearOperator<Scalar>  right_; 
    VectorSpace<Scalar> domain_;
    VectorSpace<Scalar> range_;
    bool subtraction_;

    /*  returns the next pair from the given arrays  */
    void getNext(int& k, int& LR, double& vLR,
                 Teuchos::Array<int>& ind,
                 Teuchos::Array<double>& val) const
    {
      if (k >= ind.size())
        {
          LR = -1;
        }
      else
        {
          LR = ind[k];
          vLR = val[k];
        }
      k++;
      return;
    }



    /* does bubble sort from the bottom up  */
    void sort(Teuchos::Array<int>& ind, 
              Teuchos::Array<double>& val) const
    {
      bool flip = false;
      for (int j = 0; j < ind.size();  j++)
        {
          flip = false;
          for (int i = ind.size()-1; i > j; i--)
            {
              if (ind[i] < ind[i-1])
                {
                  flip = true;
                  int t = ind[i-1];
                  ind[i-1] = ind[i];
                  ind[i] = t;

                  double v = val[i-1];
                  val[i-1] = val[i];
                  val[i] = v;
                }
            }
          if (!flip) return;
        }
    }



  };
}

#endif
