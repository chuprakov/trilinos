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

#ifndef TSFNONLINEAROPERATORBASE_HPP
#define TSFNONLINEAROPERATORBASE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFHandleable.hpp"
#include "TSFVector.hpp"

namespace TSFExtended
{
  using TSFCore::Index;
  using namespace Teuchos;

  /** 
   * Base class for nonlinear operators
   */
  template <class Scalar>
  class NonlinearOperatorBase 
    : public Handleable<NonlinearOperatorBase<Scalar> >
    {
    public:
      /** */
      NonlinearOperatorBase(const VectorSpace<Scalar>& domain,
                            const VectorSpace<Scalar>& range) 
        : domain_(domain.ptr()), range_(range.ptr()) 
      {;}
                            
      /** */
      const RefCountPtr<const TSFCore::VectorSpace<Scalar> >& domain() const 
      {return domain_;}

      /** */
      const RefCountPtr<const TSFCore::VectorSpace<Scalar> >& range() const 
      {return range_;}

      /** */
      virtual void apply(const Vector<Scalar>& in,
                         Vector<Scalar>& out) const = 0 ;

      /** */
      virtual RefCountPtr<TSFCore::LinearOp<Scalar> > jacobian(const Vector<Scalar>& x) const = 0 ;

    private:
      /** */
      RefCountPtr<const TSFCore::VectorSpace<Scalar> > domain_;

      /** */
      RefCountPtr<const TSFCore::VectorSpace<Scalar> > range_;
    };



 
}


#endif
