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

#ifndef TSFOPDESCRIBABLEBYTYPEID_HPP
#define TSFOPDESCRIBABLEBYTYPEID_HPP

#include "TSFConfigDefs.hpp"
#include "TSFDescribableByTypeID.hpp"
#include "TSFCoreLinearOp.hpp"

namespace TSFExtended
{
  /**
   * Class that extends DescribableByTypeID that is specific for
   * operators.  It gives the type and dimenstions of the operator.
   *
   * @author Paul T. Boggs (ptboggs@sandia.gov)
   */

  template <class Scalar> 
  class OpDescribableByTypeID : public DescribableByTypeID,
				public TSFCore::LinearOp<Scalar>
  {
  public:
    /** Virtual dtor */
    virtual ~OpDescribableByTypeID(){;}

    /**
     * Overwrite the describe(int depth) to append the dimensions to
     * the string.  
     *
     *@param depth int giving the number of 3 space tabs
     * at beginning of line
     */

    virtual string describe(int depth) const
    {
      string ret = "";
      for (int i = 0; i < depth; i++)
	{
	  ret.append("   ");
	}
      ret.append(typeName());
      /* The Teuchos namespace qualifier for toString()
       * seems to be required on the Intel compiler */
      ret.append(" of dimension " + Teuchos::toString(this->range()->dim()) + " by "
		 + Teuchos::toString(this->domain()->dim()));
      return ret;
    }
    

  private:

  };
}

#endif

