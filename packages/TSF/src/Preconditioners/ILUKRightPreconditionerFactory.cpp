#include "ILUKRightPreconditionerFactory.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFUtils.h"
#include "TSFLinearSolverBase.h"

using namespace TSF;
using std::string;

ILUKRightPreconditionerFactory::ILUKRightPreconditionerFactory(int fillLevels,
                                                               int overlapFill)
	: TSFPreconditionerFactoryBase(), fillLevels_(fillLevels),
		overlapFill_(overlapFill)
{}

TSFPreconditioner ILUKRightPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& A) const
{
	TSFPreconditioner rtn;

	if (A.isMatrixOperator())
		{
			const TSFSmartPtr<const TSFMatrixOperator> matrix = A.getMatrix();
			matrix->getILUKRightPreconditioner(fillLevels_, overlapFill_, rtn);
		}
	else
		{
			TSFError::raise("ILUKRightPreconditionerFactory::createPreconditioner called "
											"for a matrix-free operator");
		}
	return rtn;
}

TSFPreconditioner ILUKRightPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& S) const
{
  TSFError::raise("ILUKRightPreconditionerFactory::createPreconditioner called "
                  "for a TSFOperatorSource. Not implemented yet.");
}

string ILUKRightPreconditionerFactory::toString() const 
{
	return "ILUK factory (fill level=" + TSFUtils::toString(fillLevels_)
		+ ")";
}
