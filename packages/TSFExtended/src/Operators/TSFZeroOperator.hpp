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

#ifndef TSFZEROOPERATOR_HPP
#define TSFZEROOPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreTypes.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFOpDescribableByTypeID.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace TSFExtended
{
using namespace Teuchos;

  /** 
   * ZeroOperator is the zero operator, which maps any vector in the
   * domain space to the zero vector in the range space.
   */
  template <class Scalar> 
  class ZeroOperator : public OpDescribableByTypeID<Scalar>
  {
  public:
    /**
     * Construct the domain and range spaces.
     */
    ZeroOperator(const TSFCore::VectorSpace<Scalar>& domain, 
		 const TSFCore::VectorSpace<Scalar>& range)
      : domain_(domain.ptr()), range_(range.ptr()) {;}

    /** Virtual dtor */
    virtual ~ZeroOperator(){;}

    /** 
     * apply returns a zero vector in the range space
     */
    virtual void apply(
                       const TSFCore::ETransp            M_trans
                       ,const TSFCore::Vector<Scalar>    &x
                       ,TSFCore::Vector<Scalar>          *y
                       ,const Scalar            alpha = 1.0
                       ,const Scalar            beta  = 0.0
                       ) const 
    {
      if (beta == 0.0)
	{
	  assign(y, 0.0);
	}
      else
	{
	  Vt_S(y, beta);
	}
    }

    /** Return the domain of the operator. */
    virtual RefCountPtr< const TSFCore::VectorSpace<Scalar> > domain() const {return domain_;}
    

    /** Return the range of the operator. */
    virtual RefCountPtr< const TSFCore::VectorSpace<Scalar> > range() const {return range_;}

  protected:
    /**
     * The vector space for the range of the operator.
     */
    RefCountPtr<const TSFCore::VectorSpace<Scalar> > range_;
 
    /**
     * The vector space for the domain of the operator.
     */
    RefCountPtr<const TSFCore::VectorSpace<Scalar> > domain_;  
  };
}

#endif
