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


using namespace TSF;

#if HAVE_PETRA

AZTECSolver::AZTECSolver()
	: TSFLinearSolverBase(),
		options_(AZ_OPTIONS_SIZE),
		parameters_(AZ_PARAMS_SIZE)
{
	/* initialize the options and parameters with Aztec's defaults */
	AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));
}

AZTECSolver::AZTECSolver(const TSFHashtable<int, int>& inputOptions,
												 const TSFHashtable<int, double>& inputParameters)
	: TSFLinearSolverBase(),
		options_(AZ_OPTIONS_SIZE),
		parameters_(AZ_PARAMS_SIZE)
{
	TSFHashtable<int, int> userOptions = inputOptions;
	TSFHashtable<int, double> userParameters = inputParameters;

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

	AztecOO aztec(A, x, b);

	aztec.SetAllAztecOptions((int*) &(options_[0]));
	aztec.SetAllAztecParams((double*) &(parameters_[0]));

	int maxIters = options_[AZ_max_iter];
	double tol = parameters_[AZ_tol];

	aztec.Iterate(maxIters, tol);

	soln = xCopy;
	return true;
}

#endif



