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

#ifndef TSFLINEAROPERATORDECL_HPP
#define TSFLINEAROPERATORDECL_HPP

#include "TSFConfigDefs.hpp"
#include "TSFHandle.hpp"
#include "TSFHandleable.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVector.hpp"
 //#include "TSFVectorSpace.hpp"
#include "Teuchos_TimeMonitor.hpp"


// template <class Scalar>
//  class TransposeOperator;

// template <class Scalar>
//  class InverseOperator;

// template <class Scalar>
//  class RowAccessibleOp;


namespace TSFExtended
{
  using TSFCore::Index;
  using namespace Teuchos;

  template <class Scalar>
  class LinearSolver;

  /** 
   * User-level linear operator class
   */
  template <class Scalar>
  class LinearOperator : public Handle<TSFCore::LinearOp<Scalar> >
    {
    public:
      /** \name Constructors, Destructors, and Assignment Operators */
      //@{
      //      HANDLE_CTORS(LinearOperator<Scalar>, TSFCore::LinearOp<Scalar>);
      /** */
      LinearOperator();
      /** */
      LinearOperator(const Handleable<TSFCore::LinearOp<Scalar> >* rawPtr);
      /** */
      LinearOperator(const RefCountPtr<TSFCore::LinearOp<Scalar> >& smartPtr);
      //@}

      /** */
      virtual const VectorSpace<Scalar> domain() const ;

      /** */
      virtual const VectorSpace<Scalar> range() const ;

      /** */
      void apply(const Vector<Scalar>& in,
                 Vector<Scalar>& out) const ;

      /** */
       void applyTranspose(const Vector<Scalar>& in,
                           Vector<Scalar>& out) const ;

//       /** */
//       LinearOperator<Scalar> transpose() const ;

//       /** */
//       LinearOperator<Scalar> form() const ;
      
      
      /** Get a stopwatch for timing vector operations */
      RefCountPtr<Time>& opTimer();

      /**
       * Return a TransposeOperator.
       */
      LinearOperator<Scalar> transpose() const ; 
      
      /**
       * Return an InverseOperator.
       */
      LinearOperator<Scalar> inverse(LinearSolver<Scalar> solver = LinearSolver<Scalar>()) const ;

      RefCountPtr<LoadableMatrix<Scalar> > matrix();

      /** */     
      void getRow(const int& row, 
		  Teuchos::Array<int>& indices, 
		  Teuchos::Array<Scalar>& values) const ;
    private:
  };
}


#endif
