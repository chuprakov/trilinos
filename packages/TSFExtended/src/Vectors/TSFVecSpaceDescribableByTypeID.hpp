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

#ifndef TSFVECSPACEDESCRIBABLEBYTYPEID_HPP
#define TSFVECSPACEDESCRIBABLEBYTYPEID_HPP

#include "TSFConfigDefs.hpp"
#include "TSFDescribableByTypeID.hpp"
#include "Thyra_VectorSpaceBase.hpp"

namespace TSFExtended
{
  /**
   * This cannot be used with a vector space that already derives from
   * a Thyra::VectorSpaceBase since it then have the diamond pattern of
   * multiple inheritance
   *
   * Class that extends DescribableByTypeID that is specific for
   * VectorSpaces.  It gives the type and dimenstions of the operator.
   *
   * @author Paul T. Boggs (ptboggs@sandia.gov)
   */

  template <class Scalar> 
  class VecSpaceDescribableByTypeID : public DescribableByTypeID,
                                      public Thyra::VectorSpaceBase<Scalar>
  {
  public:
    /** Virtual dtor */
    virtual ~VecSpaceDescribableByTypeID(){;}

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
      ret.append(" of dimension " + Teuchos::Utils::toString(this->dim()));
      return ret;
    }
    

  private:

  };
}

#endif

