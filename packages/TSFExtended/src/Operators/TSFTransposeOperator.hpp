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

#ifndef TSFTRANSPOSEOPERATOR_HPP
#define TSFTRANSPOSEOPERATOR_HPP

#include "TSFConfigDefs.hpp"
 #include "TSFCoreVectorSpace.hpp"
#include "TSFOpDescribableByTypeID.hpp"
#include "TSFLinearOperatorDecl.hpp"
//#include "TSFLoadableMatrix.hpp"
//#include "TSFExplicitlyTransposeableOp.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace TSFExtended
{
  /** 
   * Class TransposeOperator represents the unformed transpose \f$A^T\f$
   * of another operator \f$A\f$. The transpose itself is never formed,
   * rather, any operator applications are carried out with the 
   * ETransp flag set to TRANS in the apply() method. 
   *
   * Should one ever want to form an explicit transpose, the 
   * form() method should be called.
   * 
   */
  template <class Scalar>
  class TransposeOperator : public OpDescribableByTypeID<Scalar>, 
			    public Handleable<TSFCore::LinearOp<Scalar> >
    //public LoadableMatrix<Scalar>
    //public ExplicitlyTransposeableOp<Scalar>
    //public Formable<Scalar>,
  {
   
  public:
    GET_RCP(TSFCore::LinearOp<Scalar>);
//     virtual RefCountPtr<TransposeOperator<Scalar> > getRcp() 
//     {return rcp(this);}





    /** Create with an operator. */
    TransposeOperator(const LinearOperator<Scalar>& op) : op_(op) {;}

    /** Virtual dtor */
    virtual ~TransposeOperator(){;}

    /** Return the domain, which is the range of the operator being
     * transposed. */
    Teuchos::RefCountPtr<const TSFCore::VectorSpace<Scalar> >  domain() const 
    {
      VectorSpace<Scalar> vs = op_.range();
      return vs.ptr();
      //return op_.range();
    }

    /** Return the range, which is the domain of the operator being
     * transposed. */
    Teuchos::RefCountPtr<const TSFCore::VectorSpace<Scalar> >  range() const 
    {
      VectorSpace<Scalar> vs = op_.domain();
      return vs.ptr();
      //return op_.domain();
    }

    /** Apply the transpose of the underlying operator */
    void apply(
               const TSFCore::ETransp            M_trans
               ,const TSFCore::Vector<Scalar>    &x
               ,TSFCore::Vector<Scalar>          *y
               ,const Scalar            alpha = 1.0
               ,const Scalar            beta = 0.0
               ) const 
    {
      op_.ptr()->apply(not_trans(M_trans), x, y, alpha, beta);
    }

    /** Form an explicit representation if supported by the underlying
     * operator. */
//     LinearOperator<Scalar> form() const 
//     {
//       const ExplicitlyTransposeable<Scalar>* et = 
//         dynamic_cast<const ExplicitlyTransposeable<Scalar>*>(op_.ptr().get());


//       TEST_FOR_EXCEPTION(et==0, runtime_error,
//                          "TransposeOperator<Scalar>::form() called where "
//                          << "the operator is unable to form an explicit "
//                          << "transpose. The operator is " 
// 			 << op_.describe());

//       return et->formTranspose();
//     }

    /** Form the transpose of this operator, which is a copy of the
     * underlying operator. */
//     LinearOperator<Scalar> formTranspose() const
//     {
//       return op_.clone();
//     }
    
  private:    
    LinearOperator<Scalar> op_;
  };
}
#endif
