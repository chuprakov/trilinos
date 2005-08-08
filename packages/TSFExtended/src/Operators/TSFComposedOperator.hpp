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

#ifndef TSFCOMPOSEDOPERATOR_HPP
#define TSFCOMPOSEDOPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "TSFSingleScalarTypeOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "TSFRowAccessibleOp.hpp"
#include "TSFHandleable.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFLinearOperatorDecl.hpp"


namespace TSFExtended
{
  /**
   * ComposedOperator is a composition of two linear operators.
   */
  template <class Scalar> 
  class ComposedOperator : public SingleScalarTypeOp<Scalar>,
                           public Handleable<SingleScalarTypeOp<Scalar> >,
                           public RowAccessibleOp<Scalar>
  {
  public:
    GET_RCP(SingleScalarTypeOp<Scalar>);
    /** 
     * Construct from a pair of linear operators.  Note: need to fix
     * the test_for_exception
     */
    ComposedOperator(const LinearOperator<Scalar>& left, 
                     const LinearOperator<Scalar>& right)
      : left_(left), 
        right_(right),
        domain_(right.domain()),
        range_(left.range())
    {
      TEST_FOR_EXCEPTION(left_.domain() != right_.range(), runtime_error,
                         "Operators used to construct CompsedOperator are "
                         << " not compatible");
    }

    /** Virtual dtor */
    virtual ~ComposedOperator(){;}

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
      Vector<Scalar> v;
      if (M_trans == Thyra::NOTRANS)
        {
          v = right_.range().createMember();
          v.zero();
          right_.generalApply(M_trans, x, v.ptr().get(), 1.0, 0.0);
          left_.generalApply(M_trans, *(v.ptr()), y, alpha, beta);
        }
      else
        {
          v = left_.range().createMember();
          v.zero();
          left_.generalApply(M_trans, x, v.ptr().get(), 1.0, 0.0);
          right_.generalApply(M_trans, *(v.ptr()), y, alpha, beta);
        }
    }

    /** Return the domain of the operator */
    virtual RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {return right_.domain().ptr();}

    /** Return the range of the operator */
    virtual RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > range() const 
    {return left_.range().ptr();}



    /** Return the kth row  */
    void getRow(const int& row, 
                Teuchos::Array<int>& indices, 
                Teuchos::Array<Scalar>& values) const
    {
      /*  The algorithm is as follows: get the "row" row of left_ and
          then compute the linear combination of rows of right_
          corresponding to the column indices of left_(row) with the
          multipliers given by the values of the row elements.  */
      indices.resize(0);
      values.resize(0);
      Teuchos::Array<int> Lind;
      Teuchos::Array<double> Lval;

      left_.getRow(row, Lind, Lval);
      sort(Lind, Lval);
      if (Lind.size() == 0)
        {
          return;
        }
      Teuchos::Array<int> Rind;
      Teuchos::Array<double> Rval;
      Rind.resize(0);
      Rval.resize(0);

      /*  The Lind[0] row of right_ isput into indices and values and
          then the values are multiplied by Lval[0].  This is the basis
          for the indexed axpy (indAxpy) to compute the rest of the
          linear combination.  */
      right_.getRow(Lind[0], indices, values);
      sort(indices, values);
      for (unsigned int i = 0; i < values.size(); i++)
        {
          values[i] = values[i] * Lval[0];
        }

      /*  Main loop to do the linear combination. */
      for (unsigned int j = 1; j < Lind.size(); j++)
        {
          Rind.resize(0);
          Rval.resize(0);
          right_.getRow(Lind[j], Rind, Rval);
          indAxpy(indices, values, Lval[j], Rind, Rval);
        }
    }





  private:

    LinearOperator<Scalar> left_;  
    LinearOperator<Scalar> right_; 
    VectorSpace<Scalar> domain_;
    VectorSpace<Scalar> range_;


    void indAxpy(Teuchos::Array<int>& indices, 
                 Teuchos::Array<double>& values, double& a,
                 Teuchos::Array<int>& Rind, 
                 Teuchos::Array<double>& Rval) const
    {
      /*  Performs an "indexed axpy."  That is it computes y = ax + y
          where y is defined by indices and values and x is defined by
          Rind and Rval.  To do this, it copies indices and values into
          Tind and Tval and rebuilds indices and values. */

      sort(indices, values);
      Teuchos::Array<int> Tind;
      Teuchos::Array<double> Tval;
      for (unsigned int i = 0; i < indices.size(); i++)
        {
          Tind.append(indices[i]);
          Tval.append(values[i]);
        }
      indices.resize(0);
      values.resize(0);
      sort(Rind, Rval);

      int kR = 0;
      int kT = 0;
      int T = 0;
      int R = 0;
      double vT = 0.0;
      double vR = 0.0;

      /* getNext gets the next pair from the specified arrays.  It
         returns R (or T) = -1 if there are no more elements. */

      getNext(kR, R, vR, Rind, Rval);
      getNext(kT, T, vT, Tind, Tval);

      /*  Main loop to rebuild indices and values. */
      int large = Rind.size() + Tind.size();
      for (int i = 0; i < large; i++)
        {
          if (R == -1 && T == -1)
            {
              break;
            }
          else if (T == -1 || (R != -1 && T > R))
            {
              indices.append(R);
              values.append(vR * a);
              getNext(kR, R, vR, Rind, Rval);
            }
          else if (R == -1 || T < R)
            {
              indices.append(T);
              values.append(vT);
              getNext(kT, T, vT, Tind, Tval);
            }
          else if (T == R)
            {
              indices.append(T);
              values.append(vT + a * vR);
              getNext(kR, R, vR, Rind, Rval);
              getNext(kT, T, vT, Tind, Tval);
            }
        }
    }




    void getNext(int& k, int& LR, double& vLR,
                 Teuchos::Array<int>& ind,
                 Teuchos::Array<double>& val) const
    {
      /* gets next pair from the given arrays.  Returns -1 in LR if no
         more elements.  */
      if ((unsigned int) k >= ind.size())
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



    void sort(Teuchos::Array<int>& ind, Teuchos::Array<double>& val) const
    {
      /* Standard bubble sort from bottom, since that's where any changes
         are made.  */
      bool flip = false;
      int len = ind.size();
      for (int j = 0; j < len;  j++)
        {
          flip = false;
          for (int i = len-1; i > j; i--)
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
