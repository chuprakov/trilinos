#include "TSFAztecSolver.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFEpetraMatrix.hpp"

#ifdef HAVE_ML
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_operator.h"
#include "ml_aztec_utils.h"
using namespace ML_Epetra;
#endif

using namespace TSFExtended;
using namespace Teuchos;





AztecSolver::AztecSolver(const std::map<int, int>& aztecOptions,
                         const std::map<int, double>& aztecParameters)
	: LinearSolverBase<double>(ParameterList()),
		options_(AZ_OPTIONS_SIZE),
		parameters_(AZ_PARAMS_SIZE),
    useML_(false),
    mlLevels_(0),
    mlSymmetric_(false),
    mlUseDamping_(false),
    mlDamping_(0.0),
    prec_(),
    aztec_recursive_iterate_(false),
    aztec_status(AZ_STATUS_SIZE),
    aztec_proc_config(AZ_PROC_SIZE)
{
  if (aztecOptions.find(AZ_ml) != aztecOptions.end()) 
    {
      useML_ = true;
    }

  if (aztecOptions.find(AZ_recursive_iterate) != aztecOptions.end()) 
    {
      aztec_recursive_iterate_ = true;
    }

#ifndef HAVE_ML
  TEST_FOR_EXCEPTION(useML_==true, runtime_error,
                     "ML is not supported in this build of TSF. "
                     "To use ML, reconfigure Trilinos with --enable-ml");
#endif

  if (aztecOptions.find(AZ_ml_levels) != aztecOptions.end()) 
    {
      mlLevels_ = aztecOptions.find(AZ_ml_levels)->second;
    }
  
  if (aztecOptions.find(AZ_ml_sym) != aztecOptions.end()) 
    {
      mlSymmetric_ = true;
    }

  if (aztecOptions.find(AZ_ml_damping) != aztecOptions.end()) 
    {
      mlUseDamping_ = true;
      mlDamping_ = aztecParameters.find(AZ_ml_damping)->second;
    }

  
	/* initialize the options and parameters with Aztec's defaults */
	AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));

	/* set user-specified options  */
  map<int, int>::const_iterator opIter;
  for (opIter=aztecOptions.begin(); opIter!=aztecOptions.end(); opIter++)
		{
      int opKey = opIter->first;
      int opValue = opIter->second;
			options_[opKey] = opValue;
		}
	
	/* set user-specified params  */
  map<int, double>::const_iterator parIter;
  for (parIter=aztecParameters.begin(); parIter!=aztecParameters.end(); 
       parIter++)
		{
      int parKey = parIter->first;
      double parValue = parIter->second;
			parameters_[parKey] = parValue;
		}
	
}




void AztecSolver::setupML(Epetra_RowMatrix* F) const
{
#ifdef HAVE_ML
  ML* ml_handle;
  ML_Aggregate* agg_object;

  ML_Set_PrintLevel(10);
  ML_Create(&ml_handle, mlLevels_);
  
  EpetraMatrix2MLMatrix(ml_handle, 0, F);

  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxCoarseSize(agg_object,30);
  ML_Aggregate_Set_Threshold(agg_object,0.0);

  if (mlSymmetric_ != true) 
    {
      ML_Aggregate_Set_DampingFactor(agg_object,0.);
    }
  mlLevels_ = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
                                                 ML_INCREASING, agg_object);
  if (mlSymmetric_ != true) 
    {
#ifdef EPETRA_MPI
      AZ_set_proc_config(&(aztec_proc_config[0]), MPI_COMM_WORLD);
#else
      AZ_set_proc_config(&(aztec_proc_config[0]), AZ_NOT_MPI);
#endif 
      
      ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, 1, 1);
    }
  else 
    {
      ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, 1, 1);
    }
  
  ML_Gen_Solver    (ml_handle, ML_MGV, 0, mlLevels_-1);
  
  MultiLevelOperator  *MLop = new MultiLevelOperator(ml_handle,
                                                     (F->OperatorDomainMap().Comm()),
                                                     (F->OperatorDomainMap()),
                                                     (F->OperatorDomainMap()));
  MLop->SetOwnership(true);
  ML_Aggregate_Destroy(&agg_object);

  prec_ = RefCountPtr<Epetra_Operator>(MLop, true);
#endif
}

SolverState<double> AztecSolver::solve(const LinearOperator<double>& op, 
                                       const Vector<double>& rhs, 
                                       Vector<double>& soln) const
{
	Vector<double> bCopy = rhs.copy();
	Vector<double> xCopy = rhs.copy();

  Epetra_Vector& b = EpetraVector::getConcrete(bCopy);
  Epetra_Vector& x = EpetraVector::getConcrete(xCopy);

	Epetra_CrsMatrix& A = EpetraMatrix::getConcrete(op);

  if (useML_) setupML(&A);

  AztecOO aztec(&A, &x, &b);
  
  aztec.SetAllAztecOptions((int*) &(options_[0]));
  aztec.SetAllAztecParams((double*) &(parameters_[0]));
  
  int maxIters = options_[AZ_max_iter];
  double tol = parameters_[AZ_tol];
  
  /* VEH/RST - check if we are using an Epetra_Operator as 
   * a preconditioner. Note that in this case user must set
   * the parameter aztec_recursive_iterate to true */
  if (prec_.get() != 0)
    aztec.SetPrecOperator(prec_.get());
  
  /* VEH/RST Parameter to check if we are calling aztec recursively.
   * If so, need to set parameter aztec_recursive_iterate to true. */
  if (aztec_recursive_iterate_)
    aztec.recursiveIterate(maxIters, tol);
  else
    aztec.Iterate(maxIters, tol);
  
  
  soln = xCopy;

  const double* status = aztec.GetAztecStatus();
  SolverStatusCode state;

  string msg;
  switch((int) status[AZ_why])
    {
    case AZ_normal:
      state = SolveConverged;
      msg = "converged";
      break;
    case AZ_param:
      state = SolveCrashed;
      msg = "failed: parameter not available";
      break;
    case AZ_breakdown:
      state = SolveCrashed;
      msg = "failed: numerical breakdown";
      break;
    case AZ_loss:
      state = SolveCrashed;
      msg = "failed: numerical loss of precision";
      break;
    case AZ_ill_cond:
      state = SolveCrashed;
      msg = "failed: ill-conditioned Hessenberg matrix in GMRES";
      break;
    case AZ_maxits:
      state = SolveFailedToConverge;
      msg = "failed: maxiters reached without converged";
      break;
    }
  SolverState<double> rtn(state, "Aztec solver " + msg, (int) status[AZ_its],
              status[AZ_r]);
  return rtn;
}





