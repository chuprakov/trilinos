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


AztecSolver::AztecSolver(const ParameterList& params)
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
  initParamMap();

  /* initialize the options and parameters with Aztec's defaults */
	AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));

  /* Set options according to the parameter list */
  ParameterList::ConstIterator iter;
  for (iter=params.begin(); iter != params.end(); ++iter)
    {
      const string& name = params.name(iter);
      const ParameterEntry& entry = params.entry(iter);
      //   cerr << "Found parameter " << name << " = " << entry << endl;
      /* Check that the param name appears in the table of Aztec params */
      if (paramMap().find(name) == paramMap().end()) continue;

      /* find the integer ID used by Aztec to identify this parameter */
      int aztecCode = paramMap()[name];

      /* We now need to figure out what to do with the value of the
       * parameter. If it is a string, then it corresponds to a
       * predefined Aztec option value. If it is an integer, then
       * it is the numerical setting for an Aztec option. If it is
       * a double, then it is the numerical setting for an Aztec
       * parameter. */
      if (entry.isType<string>())
        {
          string val = getValue<string>(entry);
          TEST_FOR_EXCEPTION(paramMap().find(val) == paramMap().end(),
                             runtime_error,
                             "Aztec solver ctor: [" << val << "] is not a "
                             "valid Aztec option value");
          int optionVal = paramMap()[val];
          options_[aztecCode] = optionVal;
        }
      else if (entry.isType<int>())
        {
          int val = getValue<int>(entry);
          options_[aztecCode] = val;
        }
      else if (entry.isType<double>())
        {
          double val = getValue<double>(entry);
          parameters_[aztecCode] = val;
        }
    }
}


AztecSolver::AztecSolver(const Teuchos::map<int, int>& aztecOptions,
                         const Teuchos::map<int, double>& aztecParameters)
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

  if (useML_) 
    {
      setupML(&A);
    }


  AztecOO aztec(&A, &x, &b);


  aztec.SetAllAztecOptions((int*) &(options_[0]));
  aztec.SetAllAztecParams((double*) &(parameters_[0]));

  aztec.CheckInput();
  
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


void AztecSolver::initParamMap()
{
  static bool first = true;
  if (first)
    {
      paramMap()["Method"]=AZ_solver;
      paramMap()["CG"]=AZ_cg;
      paramMap()["GMRES"]=AZ_gmres;
      paramMap()["CGS"]=AZ_cgs;
      paramMap()["TFQMR"]=AZ_tfqmr;
      paramMap()["BICGSTAB"]=AZ_bicgstab;
      paramMap()["Direct"]=AZ_lu;
      paramMap()["Precond"]=AZ_precond;
      paramMap()["None"]=AZ_none;
      paramMap()["Jacobi"]=AZ_Jacobi;
      paramMap()["Neumann Series"]=AZ_Neumann;
      paramMap()["Symmetric Gauss-Seidel"]=AZ_sym_GS;
      paramMap()["Least-Squares Polynomial"]=AZ_ls;
      paramMap()["Domain Decomposition"]=AZ_dom_decomp;
      paramMap()["Subdomain Solver"]=AZ_subdomain_solve;
      paramMap()["Approximate Sparse LU"]=AZ_lu;
      paramMap()["Saad ILUT"]=AZ_ilut;
      paramMap()["ILU"]=AZ_ilu;
      paramMap()["RILU"]=AZ_rilu;
      paramMap()["Block ILU"]=AZ_bilu;
      paramMap()["Incomplete Cholesky"]=AZ_icc;
      paramMap()["Residual Scaling"]=AZ_conv;
      paramMap()["Initial"]=AZ_r0;
      paramMap()["RHS"]=AZ_rhs;
      paramMap()["Matrix"]=AZ_Anorm;
      paramMap()["Solution"]=AZ_sol;
      paramMap()["No Scaling"]=AZ_noscaled;
      paramMap()["Verbosity"]=AZ_output;
      paramMap()["All"]=AZ_all;
      paramMap()["Silent"]=AZ_none;
      paramMap()["Warnings"]=AZ_warnings;
      paramMap()["Final Residual"]=AZ_last;
      paramMap()["Graph Fill"]=AZ_graph_fill;
      paramMap()["Max Iterations"]=AZ_max_iter;
      paramMap()["Polynomial Order"]=AZ_poly_ord;
      paramMap()["Overlap"]=AZ_overlap;
      paramMap()["Overlap Type"]=AZ_type_overlap;
      paramMap()["Standard"]=AZ_standard;
      paramMap()["Symmetric"]=AZ_symmetric;
      paramMap()["Restart Size"]=AZ_kspace;
      paramMap()["Reorder ILU"]=AZ_reorder;
      paramMap()["Keep Factorization"]=AZ_keep_info;
      paramMap()["GMRES Orthogonalization"]=AZ_orthog;
      paramMap()["Classical Gram-Schmidt"]=AZ_classic;
      paramMap()["Modified Gram-Schmidt"]=AZ_modified;
      paramMap()["Auxiliary Vector"]=AZ_aux_vec;
      paramMap()["Residual"]=AZ_resid;
      paramMap()["Random"]=AZ_rand;
      paramMap()["Tolerance"]=AZ_tol;
      paramMap()["Drop Tolerance"]=AZ_drop;
      paramMap()["Fill Ratio"]=AZ_ilut_fill;
      paramMap()["Damping"]=AZ_omega;

      first = false;
    }
}


