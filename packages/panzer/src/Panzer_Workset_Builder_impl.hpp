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

#ifndef PANZER_WORKSET_BUILDER_IMPL_HPP
#define PANZER_WORKSET_BUILDER_IMPL_HPP

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "Panzer_Workset.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Shards_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

// Intrepid
#include "Shards_CellTopology.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_Basis.hpp"

template<typename ArrayT>
Teuchos::RCP< std::vector<panzer::Workset> > 
panzer::buildWorksets(const panzer::PhysicsBlock& pb,
		      const std::vector<std::size_t>& local_cell_ids,
		      const ArrayT& vertex_coordinates)
{
  using std::vector;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::IntrepidFieldContainerFactory arrayFactory;

  std::size_t total_num_cells = local_cell_ids.size();

  std::size_t workset_size = pb.cellData().numCells();

  Teuchos::RCP< std::vector<panzer::Workset> > worksets_ptr = 
    Teuchos::rcp(new std::vector<panzer::Workset>);
  std::vector<panzer::Workset>& worksets = *worksets_ptr;
   
  // special case for 0 elements!
  if(total_num_cells==0) {

     // Setup integration rules and basis
     RCP<vector<int> > ir_degrees = rcp(new vector<int>(0));
     RCP<vector<string> > basis_names = rcp(new vector<string>(0));
      
     worksets.resize(1);
     std::vector<panzer::Workset>::iterator i = worksets.begin();
     i->num_cells = 0;
     i->block_id = pb.elementBlockID();
     i->ir_degrees = ir_degrees;
     i->basis_names = basis_names;

     const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = pb.getIntegrationRules();
     
     for (std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
	  ir_itr != int_rules.end(); ++ir_itr) {
         
       RCP<panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > > iv = 
	 rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
	 
       iv->setupArrays(ir_itr->second);

       ir_degrees->push_back(ir_itr->first);
       i->int_rules.push_back(iv);
     }

     const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases= pb.getBases();

     // Need to create all combinations of basis/ir pairings 
     for (std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
	  ir_itr != int_rules.end(); ++ir_itr) {

       for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
	    b_itr != bases.end(); ++b_itr) {

	 RCP<panzer::BasisIRLayout> b_layout = rcp(new panzer::BasisIRLayout(b_itr->second,*ir_itr->second));
	 
	 RCP<panzer::BasisValues<double,Intrepid::FieldContainer<double> > > bv = 
	   rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
	 
	 bv->setupArrays(b_layout,arrayFactory);

	 basis_names->push_back(b_layout->name());
	 i->bases.push_back(bv);
       }

     }

     i->details.push_back(Teuchos::rcpFromRef(*i));

     return worksets_ptr;
  } // end special case

  {
    std::size_t num_worksets = total_num_cells / workset_size;
    bool last_set_is_full = true;
    std::size_t last_workset_size = total_num_cells % workset_size;
    if (last_workset_size != 0) {
      num_worksets += 1;
      last_set_is_full = false;
    }    

    worksets.resize(num_worksets);
    std::vector<panzer::Workset>::iterator i;
    for (i = worksets.begin(); i != worksets.end(); ++i)
      i->num_cells = workset_size;
	 
    if (!last_set_is_full) {
      worksets.back().num_cells = last_workset_size;
    }
  }

  // assign workset cell local ids
  std::vector<std::size_t>::const_iterator local_begin = local_cell_ids.begin();
  for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
    std::vector<std::size_t>::const_iterator begin_iter = local_begin;
    std::vector<std::size_t>::const_iterator end_iter = begin_iter + wkst->num_cells;
    local_begin = end_iter;
    wkst->cell_local_ids.assign(begin_iter,end_iter);
    wkst->cell_vertex_coordinates.resize(workset_size,
					 vertex_coordinates.dimension(1),
					 vertex_coordinates.dimension(2));
    wkst->block_id = pb.elementBlockID();
    wkst->subcell_dim = pb.cellData().baseCellDimension();
    wkst->subcell_index = 0;
    wkst->details.push_back(Teuchos::rcpFromRef(*wkst));
  }
  
  TEUCHOS_ASSERT(local_begin == local_cell_ids.end());

  // Copy cell vertex coordinates into local workset arrays
  std::size_t offset = 0;
  for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
    for (std::size_t cell = 0; cell < wkst->num_cells; ++cell)
      for (std::size_t vertex = 0; vertex < Teuchos::as<std::size_t>(vertex_coordinates.dimension(1)); ++ vertex)
	for (std::size_t dim = 0; dim < Teuchos::as<std::size_t>(vertex_coordinates.dimension(2)); ++ dim)
	  wkst->cell_vertex_coordinates(cell,vertex,dim) = vertex_coordinates(cell + offset,vertex,dim);

    offset += wkst->num_cells;
  }

  TEUCHOS_ASSERT(offset == Teuchos::as<std::size_t>(vertex_coordinates.dimension(0)));
  
  // Set ir and basis arrayskset
  RCP<vector<int> > ir_degrees = rcp(new vector<int>(0));
  RCP<vector<string> > basis_names = rcp(new vector<string>(0));
  for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
    wkst->ir_degrees = ir_degrees;
    wkst->basis_names = basis_names;
  }

  const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = pb.getIntegrationRules();

  for (std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
       ir_itr != int_rules.end(); ++ir_itr) {
    
    ir_degrees->push_back(ir_itr->first);

    for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
      
      RCP<panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > > iv = 
	rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
    
      iv->setupArrays(ir_itr->second);
      iv->evaluateValues(wkst->cell_vertex_coordinates);
      
      wkst->int_rules.push_back(iv);

    }
  }

  const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases= pb.getBases();
     

  // Need to create all combinations of basis/ir pairings 

  // Loop over ir
  for (std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
       ir_itr != int_rules.end(); ++ir_itr) {
    
    // Loop over basis
    for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
	 b_itr != bases.end(); ++b_itr) {
      
      RCP<panzer::BasisIRLayout> b_layout = rcp(new panzer::BasisIRLayout(b_itr->second,*ir_itr->second));
      basis_names->push_back(b_layout->name());
      
      // Loop over worksets
      for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
	
	RCP<panzer::BasisValues<double,Intrepid::FieldContainer<double> > > bv = 
	  rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
	
	bv->setupArrays(b_layout,arrayFactory);
	
	std::size_t int_degree_index = 
	  std::distance(ir_degrees->begin(), 
			std::find(ir_degrees->begin(), 
				  ir_degrees->end(), 
				  ir_itr->second->order()));
	
	bv->evaluateValues(wkst->int_rules[int_degree_index]->cub_points,
			   wkst->int_rules[int_degree_index]->jac,
			   wkst->int_rules[int_degree_index]->jac_det,
			   wkst->int_rules[int_degree_index]->jac_inv,
			   wkst->int_rules[int_degree_index]->weighted_measure,
			   wkst->cell_vertex_coordinates);

	wkst->bases.push_back(bv);

      }
      
    }
    
  }

  return worksets_ptr;
}

// ****************************************************************
// ****************************************************************

template<typename ArrayT>
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
panzer::buildBCWorkset(const panzer::PhysicsBlock & volume_pb,
		       const std::vector<std::size_t>& local_cell_ids,
		       const std::vector<std::size_t>& local_side_ids,
		       const ArrayT& vertex_coordinates)
{
  using std::vector;
  using std::map;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::IntrepidFieldContainerFactory arrayFactory;

  // key is local face index, value is workset with all elements
  // for that local face
  Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets_ptr = 
    Teuchos::rcp(new std::map<unsigned,panzer::Workset>);

  // All elements of boundary condition should go into one workset.
  // However due to design of Intrepid (requires same basis for all
  // cells), we have to separate the workset based on the local side
  // index.  Each workset for a boundary condition is associated with
  // a local side for the element
  
  // std::cout << "local_side_ids.size() = " << local_side_ids.size() 
  // 	    << std::endl;

  TEUCHOS_ASSERT(local_side_ids.size() == local_cell_ids.size());
  TEUCHOS_ASSERT(local_side_ids.size() == static_cast<std::size_t>(vertex_coordinates.dimension(0)));

  // key is local face index, value is a pair of cell index and vector of element local ids
  std::map<unsigned,std::vector<std::pair<std::size_t,std::size_t> > > element_list;
  for (std::size_t cell=0; cell < local_cell_ids.size(); ++cell)
    element_list[local_side_ids[cell]].push_back(std::make_pair(cell,local_cell_ids[cell])); 

  std::map<unsigned,panzer::Workset>& worksets = *worksets_ptr;

  // create worksets 
  std::map<unsigned,std::vector<std::pair<std::size_t,std::size_t> > >::const_iterator side;
  for (side = element_list.begin(); side != element_list.end(); ++side) {

    std::vector<std::size_t>& cell_local_ids = worksets[side->first].cell_local_ids;
    Intrepid::FieldContainer<double> & coords = worksets[side->first].cell_vertex_coordinates;
    coords.resize(side->second.size(),vertex_coordinates.dimension(1),vertex_coordinates.dimension(2));
    for (std::size_t cell = 0; cell < side->second.size(); ++cell) {
      cell_local_ids.push_back(side->second[cell].second);

      for (std::size_t vertex = 0; vertex < Teuchos::as<std::size_t>(vertex_coordinates.dimension(1)); ++ vertex)
	for (std::size_t dim = 0; dim < Teuchos::as<std::size_t>(vertex_coordinates.dimension(2)); ++ dim)
	  coords(cell,vertex,dim) = vertex_coordinates(side->second[cell].first,vertex,dim);
    }
    worksets[side->first].num_cells = worksets[side->first].cell_local_ids.size();
    worksets[side->first].block_id = volume_pb.elementBlockID();
    worksets[side->first].subcell_dim = volume_pb.cellData().baseCellDimension() - 1;
    worksets[side->first].subcell_index = side->first;
    worksets[side->first].details.push_back(Teuchos::rcpFromRef(worksets[side->first]));
  }

  // setup the integration rules and bases
  for (std::map<unsigned,panzer::Workset>::iterator wkst = worksets.begin();
       wkst != worksets.end(); ++wkst) {
      
    const panzer::CellData side_cell_data(wkst->second.num_cells,
					  static_cast<int>(wkst->first),
					  volume_pb.cellData().getCellTopology());

    RCP<panzer::PhysicsBlock> side_pb = volume_pb.copyWithCellData(side_cell_data);

    const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = side_pb->getIntegrationRules();

    wkst->second.ir_degrees = rcp(new vector<int>(0));
    wkst->second.basis_names = rcp(new vector<string>(0));

    for (std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
	 ir_itr != int_rules.end(); ++ir_itr) {

      wkst->second.ir_degrees->push_back(ir_itr->first);
      
      RCP<panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > > iv = 
	rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
    
      iv->setupArrays(ir_itr->second);
      iv->evaluateValues(wkst->second.cell_vertex_coordinates);
      
      wkst->second.int_rules.push_back(iv);

      
      // Need to create all combinations of basis/ir pairings 
      const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases = side_pb->getBases();
      
      for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
	   b_itr != bases.end(); ++b_itr) {
	
	RCP<panzer::BasisIRLayout> b_layout = rcp(new panzer::BasisIRLayout(b_itr->second,*ir_itr->second));
	wkst->second.basis_names->push_back(b_layout->name());
      
	RCP<panzer::BasisValues<double,Intrepid::FieldContainer<double> > > bv = 
	  rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
	
	bv->setupArrays(b_layout,arrayFactory);
	
	std::size_t int_degree_index = 
	  std::distance(wkst->second.ir_degrees->begin(), 
			std::find(wkst->second.ir_degrees->begin(), 
				  wkst->second.ir_degrees->end(), 
				  ir_itr->second->order()));
	
	bv->evaluateValues(wkst->second.int_rules[int_degree_index]->cub_points,
			   wkst->second.int_rules[int_degree_index]->jac,
			   wkst->second.int_rules[int_degree_index]->jac_det,
			   wkst->second.int_rules[int_degree_index]->jac_inv,
			   wkst->second.int_rules[int_degree_index]->weighted_measure,
			   wkst->second.cell_vertex_coordinates);

	wkst->second.bases.push_back(bv);
      }
    }
  }

  return worksets_ptr;
}

#endif
