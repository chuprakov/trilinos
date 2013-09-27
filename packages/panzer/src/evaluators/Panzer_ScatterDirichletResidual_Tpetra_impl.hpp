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

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_TPETRA_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_TPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template<typename Traits,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, Traits,LO,GO,NodeT>::
ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
  local_side_id_ = p.get<int>("Local Side ID");

  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.get<bool>("Check Apply BC");
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT> 
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, Traits,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  fieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);

    if (checkApplyBC_)
      this->utils.setFieldData(applyBC_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, Traits,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<LOC> tpetraContainer 
          = Teuchos::rcp_dynamic_cast<LOC>(d.getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = tpetraContainer->get_x();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Residual, Traits,LO,GO,NodeT>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<GO> GIDs;
   std::vector<LO> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f(); 

   Teuchos::ArrayRCP<double> r_array = r->get1dViewNonConst();
   Teuchos::ArrayRCP<double> dc_array = dirichletCounter_->get1dViewNonConst();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!


   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = r->getMap()->getLocalElement(GIDs[i]);

      std::vector<bool> is_owned(GIDs.size(), false);
      globalIndexer_->ownedIndices(GIDs,is_owned);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
         const std::pair<std::vector<int>,std::vector<int> > & indicePair 
               = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
         const std::vector<int> & elmtOffset = indicePair.first;
         const std::vector<int> & basisIdMap = indicePair.second;
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            if(lid<0) // not on this processor!
               continue;

            int basisId = basisIdMap[basis];

	    if (checkApplyBC_)
	      if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
		continue;

            r_array[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

            // record that you set a dirichlet condition
            dc_array[lid] = 1.0;
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************


template<typename Traits,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, Traits,LO,GO,NodeT>::
ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
  local_side_id_ = p.get<int>("Local Side ID");

  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.get<bool>("Check Apply BC");
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Tangent");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT> 
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, Traits,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  fieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);

    if (checkApplyBC_)
      this->utils.setFieldData(applyBC_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, Traits,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<LOC> tpetraContainer 
          = Teuchos::rcp_dynamic_cast<LOC>(d.getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = tpetraContainer->get_x();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent, Traits,LO,GO,NodeT>::
evaluateFields(typename Traits::EvalData workset)
{ 
   TEUCHOS_ASSERT(false);

   std::vector<GO> GIDs;
   std::vector<LO> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f(); 

   Teuchos::ArrayRCP<double> r_array = r->get1dViewNonConst();
   Teuchos::ArrayRCP<double> dc_array = dirichletCounter_->get1dViewNonConst();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!


   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = r->getMap()->getLocalElement(GIDs[i]);

      std::vector<bool> is_owned(GIDs.size(), false);
      globalIndexer_->ownedIndices(GIDs,is_owned);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
         const std::pair<std::vector<int>,std::vector<int> > & indicePair 
               = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
         const std::vector<int> & elmtOffset = indicePair.first;
         const std::vector<int> & basisIdMap = indicePair.second;
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            if(lid<0) // not on this processor!
               continue;

            int basisId = basisIdMap[basis];

	    if (checkApplyBC_)
	      if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
		continue;

            r_array[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basisId).val();

            // record that you set a dirichlet condition
            dc_array[lid] = 1.0;
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, Traits,LO,GO,NodeT>::
ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
  local_side_id_ = p.get<int>("Local Side ID");
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.get<bool>("Check Apply BC");
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT> 
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, Traits,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);

    if (checkApplyBC_)
      this->utils.setFieldData(applyBC_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
  num_eq = scatterFields_.size();
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, Traits,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract linear object container
  tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.getDataObject(globalDataKey_));

  if(tpetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.getDataObject(globalDataKey_),true)->getGhostedLOC();
    tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<LOC> tpetraContainer 
          = Teuchos::rcp_dynamic_cast<LOC>(d.getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = tpetraContainer->get_x();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}

// **********************************************************************
template<typename Traits,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian, Traits,LO,GO,NodeT>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<GO> GIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> r = tpetraContainer_->get_f(); 
   Teuchos::RCP<typename LOC::CrsMatrixType> Jac = tpetraContainer_->get_A();

   Teuchos::ArrayRCP<double> r_array = r->get1dViewNonConst();
   Teuchos::ArrayRCP<double> dc_array = dirichletCounter_->get1dViewNonConst();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 
      const std::vector<LO> & LIDs = globalIndexer_->getElementLIDs(cellLocalId); 

      std::vector<bool> is_owned(GIDs.size(), false);
      globalIndexer_->ownedIndices(GIDs,is_owned);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
         const std::pair<std::vector<int>,std::vector<int> > & indicePair 
               = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
         const std::vector<int> & elmtOffset = indicePair.first;
         const std::vector<int> & basisIdMap = indicePair.second;
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            if(lid<0) // not on this processor
               continue;

            int basisId = basisIdMap[basis];

	    if (checkApplyBC_)
	      if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
		continue;

            // zero out matrix row
            {
               std::size_t sz = Jac->getNumEntriesInLocalRow(lid);
               std::size_t numEntries = 0;
               Teuchos::Array<LO> rowIndices(sz);
               Teuchos::Array<double> rowValues(sz);

               // Jac->getLocalRowView(lid,numEntries,rowValues,rowIndices);
               Jac->getLocalRowCopy(lid,rowIndices,rowValues,numEntries);

               for(std::size_t i=0;i<numEntries;i++)
                  rowValues[i] = 0.0;

               Jac->replaceLocalValues(lid,rowIndices,rowValues);
            }
 
            GO gid = GIDs[offset];
            const ScalarT & scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);
    
            r_array[lid] = scatterField.val();
            dc_array[lid] = 1.0; // mark row as dirichlet
    
            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);
    
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);
            TEUCHOS_ASSERT(jacRow.size()==GIDs.size());
    
            Jac->replaceGlobalValues(gid, GIDs, jacRow);
         }
      }
   }
}

// **********************************************************************

#endif
