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

#ifndef FEAPP_DAKOTAPCE_HPP
#define FEAPP_DAKOTAPCE_HPP

// Dakota includes

#include "system_defs.h"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "DakotaInterface.H"
//#include "DirectApplicInterface.H"

//Trilinos includes
#include "Teuchos_RCP.hpp"

//TriKota includes
#include "FEApp_DakotaElementInterface.hpp"

namespace FEApp {

  class DakotaPCE {
  public:

    // Constructor and destructor
    
    DakotaPCE(const char* dakota_in="dakota.in",  
	      const char* dakota_out="dakota.out",
	      const char* dakota_err="dakota.err",
	      const char* dakota_restart_out="dakota_restart.out");

    ~DakotaPCE() {};
    
    MPI_Comm getAnalysisComm(); 

    Dakota::ProblemDescDB& getProblemDescDB();
    
    Model getDakotaModel();
    
    void run(Teuchos::RCP<FEApp::DakotaElementInterface> elemInterface);
    
  protected:
    
    Dakota::ParallelLibrary parallel_lib;
    Dakota::ProblemDescDB problem_db;
    Dakota::Strategy selected_strategy;

  }; // end of class DakotaPCE

} // namespace FEApp

#endif //  FEAPP_DAKOTAPCE_HPP

