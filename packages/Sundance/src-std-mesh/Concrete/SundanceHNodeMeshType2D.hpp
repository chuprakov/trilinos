/* @HEADER@ */
// ************************************************************************
//
//                              Sundance
//                 Copyright (2005) Sandia Corporation
//
// Copyright (year first published) Sandia Corporation.  Under the terms
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov),
// Sandia National Laboratories, Livermore, California, USA
//
// ************************************************************************
/* @HEADER@ */

/*
 * SundanceHNodeMeshType2D.hpp
 *
 *  Created on: Mar 10, 2010
 *      Author: benk
 */

#ifndef SUNDANCE_HNODEMESHTYPE2D_HPP_
#define SUNDANCE_HNODEMESHTYPE2D_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceHNodeMesh2D.hpp"

namespace Sundance
{
  using namespace Teuchos;

  /**
   *
   */
  class HNodeMeshType2D : public MeshTypeBase
  {
  public:
    /** Empty ctor */
	  HNodeMeshType2D(const MeshEntityOrder& order=ExodusMeshOrder)
	    : order_(order) {;}

    /** virtual dtor */
    virtual ~HNodeMeshType2D(){;}

    /** Create a mesh of the given dimension */
    virtual RCP<MeshBase> createEmptyMesh(int dim,
                                          const MPIComm& comm) const
    // this line is never used since we create directly the mesh at the Mesher
    {return rcp(new HNodeMesh2D(dim, comm, order_));}

    /** */
    string description() const {return "HNodeMeshType2D";}

    /** Return a ref count pointer to self */
    virtual RCP<MeshTypeBase> getRcp() {return rcp(this);}

  private:
    MeshEntityOrder order_;

  };
}

#endif /* SUNDANCEHNODEMESHTYPE2D_HPP_ */
