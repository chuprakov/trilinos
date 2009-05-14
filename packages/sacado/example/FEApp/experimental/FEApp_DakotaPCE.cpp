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

#include "FEApp_DakotaPCE.hpp"
#include "DakotaIterator.H"

// Dakota driver when linking in as a library 
// Assumes MPI_COMM_WORLD both for Dakota and the model evaluation
FEApp::DakotaPCE::
DakotaPCE(const char* dakota_in,  
          const char* dakota_out,
          const char* dakota_err,
          const char* dakota_restart_out)
  :  
  parallel_lib(), 
  problem_db(parallel_lib)
{
  parallel_lib.specify_outputs_restart(dakota_out, dakota_err, NULL,
                                       dakota_restart_out, 0);
  problem_db.manage_inputs(dakota_in);
  problem_db.check_input();

  // instantiate the strategy
  selected_strategy = Strategy(problem_db);
}

MPI_Comm 
FEApp::DakotaPCE::
getAnalysisComm()
{
  Model& first_model = *(problem_db.model_list().begin());
  MPI_Comm analysis_comm =
     first_model.parallel_configuration_iterator()->ea_parallel_level().server_intra_communicator();

  return analysis_comm;
}

Dakota::ProblemDescDB& 
FEApp::DakotaPCE::
getProblemDescDB()
{
  return problem_db;
}

Dakota::Model 
FEApp::DakotaPCE::
getDakotaModel()
{
  //return *(problem_db.model_list().begin());
  return problem_db.iterator_list().begin()->iterated_model();
}
  
void 
FEApp::DakotaPCE::
run(Teuchos::RCP<FEApp::DakotaElementInterface> elemInterface)
{

  Dakota::Model& first_model = *(problem_db.model_list().begin());
  Dakota::Interface& interface  = first_model.interface();

  // Pass a pointer to a Dakota::DirectApplicInterface
  interface.assign_rep(elemInterface.get(), false);

  selected_strategy.run_strategy();
}
