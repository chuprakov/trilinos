// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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
// ***********************************************************************
// @HEADER

#include "MLSolverFactory.h"

using namespace Meros;
using namespace TSF;


MLSolverFactory::MLSolverFactory(bool isSymmetric)
  : isSymmetric_(isSymmetric)
{}

MLSolverFactory::~MLSolverFactory()
{}

TSFLinearSolver MLSolverFactory::createSolver(const TSFLinearOperator& F) const
{


  // #include "ml_epetra_operator.h"
  // #include "ml_aztec_utils.h"
  // int ML_TSF_defaults(TSF::TSFLinearSolver &FSolver, 
  //		    ML_solverData *solver_data,
  //		    bool symmetric, Epetra_RowMatrix *F)


  ML *ml_handle;
  ML_Aggregate *agg_object;
  ML_solverData *solver_data;
  
  int N_levels = 10;
  ML_Set_PrintLevel(10);
  ML_Create(&ml_handle, N_levels);
  solver_data->ml = ml_handle;

  
  // get the epetra matrix from F
  // put in check to exit nicely if F isn't a PetraMatrix
  Epetra_RowMatrix F_epet = (Epectra_RowMatrix)(F.getConcrete());

  // Convert F to ML matrix
  EpetraMatrix2MLMatrix(ml_handle, 0, F);
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxCoarseSize(agg_object,30);
  if (isSymmetric_ != true) {
    ML_Set_Symmetrize(ml_handle, ML_TRUE);
    ML_Aggregate_Set_DampingFactor(agg_object,0.25);
  }
  N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
                                                 ML_INCREASING, 
                                                 agg_object);
  if (isSymmetric != true) {
    if (solver_data->aztec_status.get() == 0) {
      double *dtemp = new double[AZ_STATUS_SIZE];
      solver_data->aztec_status = TSFSmartPtr<double>(dtemp, true);
    }
    if (solver_data->aztec_proc_config.get() == 0) {
      int *itemp = new int[AZ_PROC_SIZE];
      solver_data->aztec_proc_config = TSFSmartPtr<int>(itemp, true);
    }
    
    int options[AZ_OPTIONS_SIZE];
    double params[AZ_PARAMS_SIZE];

#ifdef EPETRA_MPI
    AZ_set_proc_config(solver_data->aztec_proc_config.get(), MPI_COMM_WORLD);
#else
    AZ_set_proc_config(solver_data->aztec_proc_config.get(), AZ_NOT_MPI);
#endif

    AZ_defaults(options, params);
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_ilut;
    options[AZ_overlap] = 1;
    params[AZ_ilut_fill] = 2.;
    
    
    for (int level = 0; level < N_levels; level++) {
      ML_Gen_SmootherAztec(ml_handle, level, options, params, 
                           solver_data->aztec_proc_config.get(), 
                           solver_data->aztec_status.get(), 
                           AZ_ONLY_PRECONDITIONER,  ML_BOTH, NULL);
    }
  }
  
  
  else 
    ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, 1, 1);
  
  ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);
  
  Epetra_ML_Operator  *MLop = new Epetra_ML_Operator(ml_handle,
                                                     (F->OperatorDomainMap().Comm()),
                                                     (F->OperatorDomainMap()),
                                                     (F->OperatorDomainMap()));
  MLop->SetOwnership(true);
  ML_Aggregate_Destroy(&agg_object);
  
  if (isSymmetric_ == true)
    solver_data->azOptions.put(AZ_solver, AZ_cg);
  else
    solver_data->azOptions.put(AZ_solver, AZ_gmres);
  
  solver_data->azOptions.put(AZ_kspace, 50);
  solver_data->azOptions.put(AZ_conv, AZ_r0);
  solver_data->azParams.put(AZ_tol, 1e-8);
  solver_data->azOptions.put(AZ_max_iter, 50);
  solver_data->azOptions.put(AZ_output, 1);
  TSFSmartPtr<Epetra_Operator> Smart_MLprec = TSFSmartPtr<Epetra_Operator>(MLop, true);
  
  TSFLinearSolver FSolver = new AZTECSolver(solver_data->azOptions, 
                                            solver_data->azParams, 
                                            Smart_MLprec);
  
  return FSolver;
}
  
