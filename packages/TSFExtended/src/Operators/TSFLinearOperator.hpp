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

#ifndef TSFLINEAROPERATOR_HPP
#define TSFLINEAROPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFHandle.hpp"
#include "TSFHandleable.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFLoadableMatrix.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "TSFVector.hpp"
#include "TSFVectorSpace.hpp"

namespace TSFExtended
{
  using TSFCore::Index;
  using namespace Teuchos;

  /** 
   * User-level linear operator class
   */
  template <class Scalar>
  class LinearOperator : public Handle<TSFCore::LinearOp<Scalar> >
    {
    public:
      /** \name Constructors, Destructors, and Assignment Operators */
      //@{
      HANDLE_CTORS(LinearOperator<Scalar>, TSFCore::LinearOp<Scalar>);
      
      //@}

      /** */
      VectorSpace<Scalar> domain() const 
      {return ptr()->domain();}

      /** */
      VectorSpace<Scalar> range() const 
      {return ptr()->range();}

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
      
      
      /** Get a stopwtach for timing vector operations */
      static RefCountPtr<Time>& opTimer()
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Low-level vector operations");
        return rtn;
      }

      RefCountPtr<LoadableMatrix<Scalar> > matrix()
      {
        RefCountPtr<LoadableMatrix<Scalar> > rtn 
          = rcp_dynamic_cast<LoadableMatrix<Scalar> >(ptr());
        return rtn;
      }
    private:
    };



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
     * create a new vector in the domain space */
    if (out.ptr().get()==0)
      {
        out = domain().createMember();
      }
    ptr()->apply(TSFCore::TRANS, *(in.ptr().get()),
                 out.ptr().get());
  }
}


#endif
