#include "AZTECSolver.h"
#include "TSFUtils.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include "TSFError.h"
#include "TSFOut.h"
#include "TSFMPI.h"
#include "PetraMatrix.h"
#include "PetraVector.h"

#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_operator.h"
#include "ml_aztec_utils.h"

using namespace TSF;



AZTECSolver::AZTECSolver()
	: TSFLinearSolverBase(),
		options_(AZ_OPTIONS_SIZE),
		parameters_(AZ_PARAMS_SIZE),
    useML_(false),
    mlLevels_(0),
    prec_(0)
{
	/* initialize the options and parameters with Aztec's defaults */
	AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));
}

AZTECSolver::AZTECSolver(const TSFHashtable<int, int>& inputOptions,
												 const TSFHashtable<int, double>& inputParameters)
	: TSFLinearSolverBase(),
		options_(AZ_OPTIONS_SIZE),
		parameters_(AZ_PARAMS_SIZE),
    useML_(false),
    mlLevels_(0),
    mlSymmetric_(false),
    mlUseDamping_(false),
    mlDamping_(0.0),
    prec_(0),
    aztec_status(AZ_STATUS_SIZE),
    aztec_proc_config(AZ_PROC_SIZE)
{
	TSFHashtable<int, int> userOptions = inputOptions;
	TSFHashtable<int, double> userParameters = inputParameters;

  /* AZ_ml is not a standard AZ option so we need to remove it from the list
   * before sending options to AZ */
  if (userOptions.containsKey(AZ_ml)) 
    {
      useML_ = true;
      userOptions.remove(AZ_ml);
      TSFOut::println("AZTECSolver is using ML");
    }

  /* AZ_ml_levels is not a standard AZ option so we need to remove it from the list
   * before sending options to AZ */
  if (userOptions.containsKey(AZ_ml_levels)) 
    {
      mlLevels_ = userOptions.get(AZ_ml_levels);
      userOptions.remove(AZ_ml_levels);
      TSFOut::printf("AZTECSolver is using %d levels of ML\n", mlLevels_);
    }
  
  /* AZ_ml_sym is not a standard AZ option so we need to remove it from the list
   * before sending options to AZ */
  if (userOptions.containsKey(AZ_ml_sym)) 
    {
      mlSymmetric_ = true;
      userOptions.remove(AZ_ml_sym);
      TSFOut::println("AZTECSolver is using symmetric ML");
    }

  /* AZ_ml_damping is not a standard AZ option so we need to remove it from the list
   * before sending options to AZ */
  if (userParameters.containsKey(AZ_ml_damping)) 
    {
      mlUseDamping_ = true;
      mlDamping_ = userParameters.get(AZ_ml_damping);
      userOptions.remove(AZ_ml_levels);
      TSFOut::printf("AZTECSolver is using ML with damping=%g\n", mlDamping_);
    }

  

	/* initialize the options and parameters with Aztec's defaults */
	AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));

	/* set user-specified options  */
	TSFArray<int> optionKeys;
	TSFArray<int> optionVals;
	userOptions.arrayify(optionKeys, optionVals);
	for (int i=0; i<optionKeys.length(); i++)
		{
			options_[optionKeys[i]] = optionVals[i];
		}
	
	/* set user-specified params  */
	TSFArray<int> paramKeys;
	TSFArray<double> paramVals;
	userParameters.arrayify(paramKeys, paramVals);
	for (int i=0; i<paramKeys.length(); i++)
		{
			parameters_[paramKeys[i]] = paramVals[i];
		}

}



AZTECSolver::~AZTECSolver(){;}

void AZTECSolver::setupML(Epetra_RowMatrix* F) const
{
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
      ML_Set_Symmetrize(ml_handle, ML_TRUE);
      if (mlUseDamping_) ML_Aggregate_Set_DampingFactor(agg_object,mlDamping_);
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

      int options[AZ_OPTIONS_SIZE];
      double params[AZ_PARAMS_SIZE];
      AZ_defaults(options, params);
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilut;
      options[AZ_overlap] = 1;
      params[AZ_ilut_fill] = 2.;
      
      for (int level = 0; level < mlLevels_; level++) {
        ML_Gen_SmootherAztec(ml_handle, level, options, params, 
                             &(aztec_proc_config[0]), &(aztec_status[0]), 
                             AZ_ONLY_PRECONDITIONER,  ML_BOTH, NULL);
      }
    }
  else 
    {
      ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, 1, 1);
    }

  ML_Gen_Solver    (ml_handle, ML_MGV, 0, mlLevels_-1);

  Epetra_ML_Operator  *MLop = new Epetra_ML_Operator(ml_handle,
                                                     (F->OperatorDomainMap().Comm()),
                                                     (F->OperatorDomainMap()),
                                                     (F->OperatorDomainMap()));
  MLop->SetOwnership(true);
  ML_Aggregate_Destroy(&agg_object);

  prec_ = TSFSmartPtr<Epetra_Operator>(MLop, true);
}

bool AZTECSolver::solve(const TSFLinearOperator& op, 
												const TSFVector& rhs, 
												TSFVector& soln) const
{
	TSFVector bCopy = rhs.copy();
	TSFVector xCopy = rhs.copy();
	const PetraVector& bpv = PetraVector::getConcrete(bCopy);
	PetraVector& xpv = PetraVector::getConcrete(xCopy);

	Epetra_Vector* b = (Epetra_Vector*) &(bpv.values());
	Epetra_Vector* x = (Epetra_Vector*) &(xpv.values());

	Epetra_CrsMatrix* A = PetraMatrix::getConcrete(op);

  if (useML_) setupML(A);

	AztecOO aztec(A, x, b);

	aztec.SetAllAztecOptions((int*) &(options_[0]));
	aztec.SetAllAztecParams((double*) &(parameters_[0]));

	int maxIters = options_[AZ_max_iter];
	double tol = parameters_[AZ_tol];

  /* VEH/RST - need to use recursiveIterate if we are using */
  /* an Epetra_Operator as a preconditioner */
	if (prec_.get() != 0) {
    aztec.SetPrecOperator(prec_.get());
    aztec.recursiveIterate(maxIters, tol);
  }
  
  else
  	aztec.Iterate(maxIters, tol);
  
	soln = xCopy;
	return true;
}





