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

#ifndef TSFLINEAROPERATORIMPL_HPP
#define TSFLINEAROPERATORIMPL_HPP

#include "TSFConfigDefs.hpp"
#include "TSFLinearOperatorDecl.hpp"
 //#include "TSFHandle.hpp"
 //#include "TSFHandleable.hpp"
 //#include "TSFCoreLinearOp.hpp"
 //#include "TSFLoadableMatrix.hpp"
 //#include "TSFRowAccessibleOp.hpp"
 //#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RefCountPtr.hpp"
//#include "TSFVector.hpp"
//#include "TSFVectorSpace.hpp"
//#include "TSFTransposeOperator.hpp"
//#include "TSFInverseOperator.hpp"
#include "TSFComposedOperator.hpp"
#include "TSFIdentityOperator.hpp"
//#include "TSFLinearSolverBaseDecl.hpp"

// template <class Scalar>
//  class TransposeOperator;
#include "TSFTransposeOperator.hpp"

template <class Scalar>
 class InverseOperator;

// template <class Scalar>
//  class IdentityOperator;


template <class Scalar>
 class RowAccessibleOp;

using namespace TSFExtended;

template <class Scalar>
LinearOperator<Scalar>::LinearOperator() : Handle<TSFCore::LinearOp<Scalar> >() {;}

// template <class Scalar>
// LinearOperator<Scalar>::LinearOperator(const LinearOperator<Scalar>* ptr)
// {
//   RefCountPtr<TSFCore::LinearOp> ptr_ = ptr->ptr();
// }


template <class Scalar>
LinearOperator<Scalar>::LinearOperator(Handleable<TSFCore::LinearOp<Scalar> >* rawPtr) 
  : Handle<TSFCore::LinearOp<Scalar> >(rawPtr) {;}

template <class Scalar>
LinearOperator<Scalar>::LinearOperator(const RefCountPtr<TSFCore::LinearOp<Scalar> >& smartPtr) 
  : Handle<TSFCore::LinearOp<Scalar> >(smartPtr) {;}

template <class Scalar>
VectorSpace<Scalar> LinearOperator<Scalar>::domain() const 
{return ptr()->domain();}

template <class Scalar>
VectorSpace<Scalar> LinearOperator<Scalar>::range() const 
{return ptr()->range();}

template <class Scalar> inline 
void LinearOperator<Scalar>::apply(const Vector<Scalar>& in,
				   Vector<Scalar>& out) const
{
  /* the result vector might not be initialized. If it's null,
   * create a new vector in the range space */
  if (out.ptr().get()==0)
    {
      out = range().createMember();
    }
  ptr()->apply(TSFCore::NOTRANS, *(in.ptr().get()),
	       out.ptr().get());
}

template <class Scalar> inline 
void LinearOperator<Scalar>::applyTranspose(const Vector<Scalar>& in,
                                            Vector<Scalar>& out) const
{
  /* the result vector might not be initialized. If it's null,
   * create a new vector in the domain space (i.e., the range space
   * of the transpose operator */
  if (out.ptr().get()==0)
    {
      out = domain().createMember();
    }
  ptr()->apply(TSFCore::TRANS, *(in.ptr().get()),
	       out.ptr().get());
}

//       LinearOperator<Scalar> form() const ;

template <class Scalar>
RefCountPtr<Time>& LinearOperator<Scalar>::opTimer()
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("Low-level vector operations");
  return rtn;
}

template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::transpose() const
{
//   TSFCore::LinearOp<Scalar>* t = 
//     new TransposeOperator<Scalar>(ptr());
  
//   Teuchos::RefCountPtr<TSFCore::LinearOp<Scalar> > rt = 
//     rcp_dynamic_cast<TSFCore::LinearOp<Scalar> >(ptr());
//   //rp(t);
  
//   LinearOperator<Scalar> op(rt);
//   return op;
  LinearOperator<Scalar> op = new TransposeOperator<Scalar>(*this);
//   RefCountPtr<TSFCore::LinearOp<double> > op = 
//     rcp(new TransposeOperator<Scalar>(*this));
// 				      //(ptr()->range())));
  return op;
}

template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::inverse(const LinearSolver<Scalar>& solver) const
{
  LinearOperator<Scalar> op = new InverseOperator<Scalar>(*this, solver);
  return op;
}

template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::operator*(const LinearOperator<Scalar>& other) const
{
  LinearOperator<Scalar> op = new ComposedOperator<Scalar>(*this, other);
  return op;
}



template <class Scalar>
RefCountPtr<LoadableMatrix<Scalar> > LinearOperator<Scalar>::matrix()
{
  RefCountPtr<LoadableMatrix<Scalar> > rtn 
    = rcp_dynamic_cast<LoadableMatrix<Scalar> >(ptr());
  return rtn;
}

template <class Scalar>
void LinearOperator<Scalar>::getRow(const int& row, 
				    Teuchos::Array<int>& indices, 
				    Teuchos::Array<Scalar>& values) const
{
  RowAccessibleOp<Scalar> val = dynamic_cast<RowAccessibleOp<Scalar> >(ptr());
  TEST_FOR_EXCEPTION(val != 0, runtime_error, 
		     "LinearOperator<Scalar>::getRow() not defined for current operator.");
  ptr()->getRow(row, indices, values);
}

#endif
