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

#ifndef PANZER_FIELD_MANAGER_BUILDER_HPP
#define PANZER_FIELD_MANAGER_BUILDER_HPP

#include <iostream>
#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_WorksetContainer.hpp"

// Forward Declarations
namespace panzer {
  struct Traits;
  struct Workset;
  template <typename LO, typename GO> class ConnManager;
  template <typename LO, typename GO> class DOFManager;
  class EquationSetFactory;
  class BCStrategyFactory;
  class PhysicsBlock;
}

namespace PHX {
  template<typename T> class FieldManager;
}  

namespace panzer {

  template <typename LO, typename GO>
  class FieldManagerBuilder {

  public:

    typedef std::map<unsigned,panzer::Workset> BCFaceWorksetMap;

    void print(std::ostream& os) const;

    const 
      std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >&
      getVolumeFieldManagers() const {return phx_volume_field_managers_;}

    const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& 
      getWorksets() const {return worksets_;}

    const std::map<panzer::BC, 
		   std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
		   panzer::LessBC>& 
      getBCFieldManagers() const {return bc_field_managers_;}

    const std::map<panzer::BC,
		   Teuchos::RCP<std::map<unsigned,panzer::Workset> >,
		   panzer::LessBC>&
      getBCWorksets() const {return bc_worksets_;}

    // The intention of the next set of functions is to simplify and eventually
    // replace the setup routine above. Its not clear that these functions
    // belong in the field manager builder. Additionally this will add increased
    // flexibility to the field manager build in that the DOFManager will be
    // specified in a more flexable and generic way. Onward.... (ECC - 1/13/11)

    /** Setup the volume field managers. This uses the passed in <code>dofManager</code>
      * and sets it for permenant use.
      */
    void setupVolumeFieldManagers(WorksetContainer & wkstContainer, 
                                  const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
				  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
				  const Teuchos::ParameterList& closure_models,
                                  const LinearObjFactory<panzer::Traits> & lo_factory,
				  const Teuchos::ParameterList& user_data);

    /** Build the BC field managers.
      */
    void setupBCFieldManagers(WorksetContainer & wkstContainer,
                              const std::vector<panzer::BC> & bcs,
                              const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
	                      const panzer::EquationSetFactory & eqset_factory,
			      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                              const panzer::BCStrategyFactory& bc_factory,
			      const Teuchos::ParameterList& closure_models,
                              const LinearObjFactory<panzer::Traits> & lo_factory,
			      const Teuchos::ParameterList& user_data);

    void writeVolumeGraphvizDependencyFiles(std::string filename_prefix,
					    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const;

  private:

    //! Phalanx volume field managers for each element block.
    std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >
      phx_volume_field_managers_;
    
    //! Volume fill worksets for each element block.
    std::vector< Teuchos::RCP<std::vector<panzer::Workset> > > worksets_;

    /*! \brief Field managers for the boundary conditions

        key is a panzer::BC object.  value is a map of
        field managers where the key is the local side index used by
        intrepid
    */
    std::map<panzer::BC, 
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC> bc_field_managers_;

    /*! \brief Worksets for the boundary conditions

        key is a panzer::BC object.  value is a map of
        worksets where the key is the local side index used by
        intrepid.  All elemenst of a boundary are in one workset for
        each local side of the sideset.
    */
    std::map<panzer::BC,
      Teuchos::RCP<std::map<unsigned,panzer::Workset> >,
      panzer::LessBC> bc_worksets_;

  };

template<typename LO, typename GO>
std::ostream& operator<<(std::ostream& os, const panzer::FieldManagerBuilder<LO,GO>& rfd);

} // namespace panzer

#include "Panzer_FieldManagerBuilder_impl.hpp"

#endif
