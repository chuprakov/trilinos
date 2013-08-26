// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

#include "Panzer_STK_SurfaceNodeNormals.hpp"

// #include "Panzer_STK_SetupUtilities.hpp"
// #include "Panzer_WorksetContainer.hpp"
// #include "Panzer_Workset_Builder.hpp"
// #include "Panzer_STK_WorksetFactory.hpp"
// #include "Panzer_CellData.hpp"

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

// #include "Shards_CellTopology.hpp"
// #include "Intrepid_FunctionSpaceTools.hpp"
// #include "Intrepid_CellTools.hpp"

namespace panzer {
  
  TEUCHOS_UNIT_TEST(node_normals, stk_testing)
  {
    using Teuchos::RCP;
    
    //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::RCP<Teuchos::FancyOStream> pout= Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    pout->setShowProcRank(true);

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("Z Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",2);
    pl->set("Z Elements",2);
    
    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);


    // This is testing for learning about stk mesh
      
    std::string sideName = "top";
    std::string blockName = "eblock-0_0_0";
    std::vector<stk::mesh::Entity*> sides;
    
    // This is for local sides, no ghosted  
    //mesh->getMySides(sideName,blockName,sides);
    RCP<stk::mesh::fem::FEMMetaData> metaData = mesh->getMetaData();
    RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();
    
    stk::mesh::Part * sidePart = mesh->getSideset(sideName);
    stk::mesh::Part * elmtPart = mesh->getElementBlockPart(blockName);
    stk::mesh::Selector side = *sidePart;
    stk::mesh::Selector block = *elmtPart;
    stk::mesh::Selector ownedBlock = metaData->universal_part() & block & side;
    //stk::mesh::Selector ownedBlock = metaData->locally_owned_part() & block & side;
    
    stk::mesh::get_selected_entities(ownedBlock,bulkData->buckets(mesh->getSideRank()),sides);
    //stk::mesh::Part* sidePart = metaData_->get_part(sideName);
    //stk::mesh::Selector side = *sidePart;
    
    std::cout << std::endl;
    
    //comm->barrier();
    
    for (std::vector<stk::mesh::Entity*>::const_iterator side=sides.begin(); side != sides.end(); ++side) {
      *pout << "side element: rank(" << (*side)->entity_rank() << ")"
	    << ", gid(" << (*side)->identifier() << ")"
	    << ", owner_rank(" << (*side)->owner_rank()
	    << ")" << std::endl;
      
      // get node relations
      std::vector<stk::mesh::Entity*> nodes;
      
      stk::mesh::PairIterRelation node_relations = (*side)->relations(mesh->getNodeRank());
      stk::mesh::PairIterRelation parent_element_relations = (*side)->relations(mesh->getElementRank());
      
      *pout << "parent element relation: " 
	    << "size(" << parent_element_relations.size() << ")"  
	    << ", entity_rank(" << parent_element_relations[0].entity_rank() << ")" 
	    << ", topo map id(" << parent_element_relations[0].identifier() << ")" 
	    << ", gid(" << parent_element_relations[0].entity()->identifier() << ")" 
	    << std::endl;
      
      pout->pushTab(4);
      for (stk::mesh::PairIterRelation::iterator node = node_relations.begin(); node != node_relations.end(); ++node) {
	*pout << "face to node relation: "
	      << "gid(" << node->entity()->identifier() << ")" 
	      << ", topo map id(" << node->identifier() << ")" 
	      << std::endl;
      }
      pout->popTab();
      
    }
    
  }
  


  TEUCHOS_UNIT_TEST(node_normals, 3D)
  {
    using Teuchos::RCP;
    
    std::cout << std::endl;

    //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::RCP<Teuchos::FancyOStream> pout= Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    pout->setShowProcRank(true);

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("Z Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",2);
    pl->set("Z Elements",2);
    
    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    std::string sideName = "top";
    std::string blockName = "eblock-0_0_0";
    
    std::map<unsigned,std::vector<double> > normals;
    
    panzer_stk::computeSidesetNodeNormals(normals,mesh,sideName,blockName,&std::cout,pout.get());

    for (std::map<unsigned,std::vector<double> >::const_iterator node = normals.begin();
	 node != normals.end(); ++node) {
      double tol = 100.0 * Teuchos::ScalarTraits<double>::eps();
      TEST_FLOATING_EQUALITY(normals[node->first][0], 0.0, tol);
      TEST_FLOATING_EQUALITY(normals[node->first][1], 1.0, tol);
      TEST_FLOATING_EQUALITY(normals[node->first][2], 0.0, tol);
    }

  }

  TEUCHOS_UNIT_TEST(node_normals, 3D_NoDebug)
  {
    using Teuchos::RCP;
    
    std::cout << std::endl;

    //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::RCP<Teuchos::FancyOStream> pout= Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    pout->setShowProcRank(true);

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("Z Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",2);
    pl->set("Z Elements",2);
    
    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    std::string sideName = "top";
    std::string blockName = "eblock-0_0_0";
    
    std::map<unsigned,std::vector<double> > normals;
    
    panzer_stk::computeSidesetNodeNormals(normals,mesh,sideName,blockName);

    for (std::map<unsigned,std::vector<double> >::const_iterator node = normals.begin();
	 node != normals.end(); ++node) {
      double tol = 100.0 * Teuchos::ScalarTraits<double>::eps();
      TEST_FLOATING_EQUALITY(normals[node->first][0], 0.0, tol);
      TEST_FLOATING_EQUALITY(normals[node->first][1], 1.0, tol);
      TEST_FLOATING_EQUALITY(normals[node->first][2], 0.0, tol);
    }

  }


}
