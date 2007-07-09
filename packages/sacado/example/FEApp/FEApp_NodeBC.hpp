// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_NODEBC_HPP
#define FEAPP_NODEBC_HPP

#include <vector>
#include "Epetra_Map.h"
#include "FEApp_AbstractNodeBCStrategy.hpp"
#include "FEApp_TemplateTypes.hpp"

namespace FEApp {

  /*!
   * \brief Concrete object to represent a nodal boundary condition
   */
  class NodeBC {
  public:

    //! Constructor
    template <typename BuilderT>
    NodeBC(const Epetra_Map& dofMap, 
	   unsigned int gid,
	   unsigned int neqn,
	   const BuilderT& builder) : global_node_id(gid) {
      strategyTM.buildObjects(builder);
      is_owned = dofMap.MyGID(gid*neqn);
    }

    //! Destructor
    virtual ~NodeBC() {}

    //! Get global node ID
    unsigned int getNodeGID() const { return global_node_id; }

    //! Get residual offsets
    const std::vector<unsigned int>& getOffsets() const { 
      return strategyTM.begin()->getOffsets(); }

    //! Return if BC is locally owned
    bool isOwned() const { return is_owned; }

    //! Get strategy
    template <typename ScalarT>
    Teuchos::RCP< FEApp::AbstractNodeBCStrategy<ScalarT> >
    getStrategy() { return strategyTM.getAsObject<ScalarT>(); }

  private:

    //! Private to prohibit copying
    NodeBC(const NodeBC&);

    //! Private to prohibit copying
    NodeBC& operator=(const NodeBC&);

  protected:

    //! Global ID of node for this BC
    unsigned int global_node_id;

    //! Template manager storing instantiations of strategies
    AbstractNodeBCStrategy_TemplateManager<ValidTypes> strategyTM;

    //! Is BC owned on this proc
    bool is_owned;

  };

}

#endif // FEAPP_NODEBC_HPP
