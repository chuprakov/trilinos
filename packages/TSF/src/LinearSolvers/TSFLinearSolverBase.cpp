#include "TSFLinearSolverBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFPreconditionerFactory.h"
#include "TSFPreconditionerFactoryBase.h"

using namespace TSF;


TSFLinearSolverBase::TSFLinearSolverBase(const TSFParameterList& params)
  : verbosity_(1), params_(defaultParameters().overrideWith(params))
{;}

TSFLinearSolverBase::~TSFLinearSolverBase()
{;}

TSFParameterList TSFLinearSolverBase::defaultParameters() const 
{
  return TSFParameterList("Empty");
}

