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

#ifndef TSFGHOSTVIEW_HPP
#define TSFGHOSTVIEW_HPP

#include "TSFAccessibleVector.hpp"
#include "TSFVector.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * GhostView is an interface for read-only views
   * of vector elements including selected
   * off-processor elements. GhostView has no standard constructor; subclasses
   * should be constructed using the importView() method of GhostImporter.
   */
  template <class Scalar>
  class GhostView : public AccessibleVector<Scalar>
  {
  public:
    /** Virtual dtor */
    virtual ~GhostView(){;}
    
    /** Indicate whether the value at the given global index is accessible
     * in this view. */
    virtual bool isAccessible(int globalIndex) const = 0 ;

  private:
  };

}

#endif
